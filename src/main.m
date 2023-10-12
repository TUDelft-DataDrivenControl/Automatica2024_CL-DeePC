%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Closed-loop Data-Enabled Predictive Control
%           Authors: R. Dinkla, S.P. Mulders, T. Oomen, J.W. van Wingerden
%           Submission to Automatica 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
yalmip('clear');
clear;
rng default;
clc;

%% Edit path - add dependencies
[mfilePath,~,~] = fileparts(mfilename('fullpath')); % get dir of this file
cd(mfilePath); % change cwd to folder containing this file

addpath(genpath("../data"),'-begin')
addpath(genpath("../config"),'-begin')

% read directories to add to path
fid = fopen('../config/dir_locs.txt');
dirs = textscan(fid,'%s','delimiter','\n');
dirs = dirs{:};
for k1=1:length(dirs)
    addpath(genpath(dirs{k1}),'-begin')
end

% add directory of this file to path
addpath(genpath(pwd),'-begin');

%% Simulation settings
model_Favoreel1999 % loads model from Favoreel 1999 - original SPC paper

% controller settings
p = 20;
f = 20;
Nmin = (p+1)*nu + p*ny; %N >= n + (p+1)*nu + p*ny
N    = 500;%Nmin;
Nbar = p+N;
Qk = 100;
Rk = 0;
dRk= 10;

% simulation length
OL_steps  = 1200; % 1200
sim_steps = 1800; % 1800
num_steps = OL_steps+sim_steps;

% predefine inputs, outputs, states, reference
c1.u = nan(nu,num_steps);
c1.y = nan(ny,num_steps);
c1.x = nan(nx,num_steps+1);
r = nan(ny,sim_steps+f-1); % +f-1 needed for simulation end
r = (-square((0:size(r,2)-1)*2*pi/(500)))*50+1*50;

% noise
Re = 0.01*eye(ny);                       % variance
e  = mvnrnd(zeros(ny,1),Re,num_steps).'; % trajectory realization

%% initial open loop simulation
x0 = zeros(nx,1); % initial state
Ru = 1*eye(nu);   % noise variance
c1.u(:,1:OL_steps)  = mvnrnd(zeros(nu,1),Ru,OL_steps).';

% simulate system
[y_ol,~,x_ol] = lsim(plant,[c1.u(:,1:OL_steps);e(:,1:OL_steps)],[],x0);
c1.y(:,1:OL_steps) = y_ol.';
c1.x(:,1:OL_steps) = x_ol.';
c1.x(:,OL_steps+1) = plant.A*c1.x(:,OL_steps) + plant.B*[c1.u(:,OL_steps); e(:,OL_steps)];

% normalization factors
u_fac = max(abs(c1.u(:,1:OL_steps)),[],2);
y_fac = max(abs(c1.y(:,1:OL_steps)),[],2);

% normalize data
c1.u = c1.u./u_fac;
c1.y = c1.y./y_fac;
r    = r./y_fac;

%
% x+ = A x + B u + K e
% y  = C x + D u + e
%
% ub = u/fu <-> u = ub*fu
% yb = y/fy <-> y = yb*fy
%
% x+ =      A x +    (B*fu) ub + K e
% yb = (C/fy) x + (D/fy*fb) ub + (1/fy) e

% normalize plant
plant.B(:,1:nu)     =             plant.B(:,1:nu)*diag(u_fac);
plant.C             = diag(y_fac)\plant.C;
plant.D(:,1:nu)     = diag(y_fac)\plant.D(:,1:nu)*diag(u_fac);
plant.D(:,nu+1:end) = diag(y_fac)\plant.D(:,nu+1:end);

%% closed-loop operation
du = 0*mvnrnd(zeros(nu,1),Ru,sim_steps).';
du_max = 3.75;
du(abs(du)>du_max) = sign(du(abs(du)>du_max))*du_max;
du = du./u_fac;

% initialize data structures for controllers
c1.uf_k    = cell(sim_steps,1);
c1.yfhat_k = c1.uf_k;
c1.G       = c1.uf_k;

% create user-defined constraints
con = struct();
con.uf = sdpvar(nu,f,'full');
con.u0 = sdpvar(nu,1,'full');
con.expr = [con.uf <=  15./u_fac;
            con.uf >= -15./u_fac;
            con.u0 - con.uf(:,1) <=  du_max./u_fac;
            con.u0 - con.uf(:,1) >= -du_max./u_fac;
            con.uf(:,1:end-1)-con.uf(:,2:end) <=  du_max./u_fac;
            con.uf(:,1:end-1)-con.uf(:,2:end) >= -du_max./u_fac];

% initialize controllers after k1 = Nbar
range_OL = OL_steps-Nbar+1:OL_steps;
Cont1 = CL_DeePC(c1.u(:,range_OL),c1.y(:,range_OL),p,f,N,Qk,Rk,dRk,use_IV=true,constr=con);

% set counters
k1 = OL_steps + 1;
k2 = 1;

% solve
[c1.uf_k{k2},c1.yfhat_k{k2},c1.G{k2}] = Cont1.solve(rf=r(:,k2:k2+f-1));
c1.u(:,k1) = c1.uf_k{k2}(:,1) + du(:,k2);
c1 = step_plant(c1,e,plant,k1);

try
    for k1 = OL_steps+2:num_steps
        k2 = k2 + 1;
        
        [c1.uf_k{k2},c1.yfhat_k{k2},c1.G{k2}] = Cont1.step(c1.u(:,k1-1),c1.y(:,k1-1), rf=r(:,k2:k2+f-1)); 
        c1.u(:,k1) = c1.uf_k{k2}(:,1) + du(:,k2);
        c1 = step_plant(c1,e,plant,k1);
    end
catch
    k2 = k2 -1;
end
%%
% reverse normalization
c1.u = c1.u.*u_fac;
c1.y = c1.y.*y_fac;
r = r.*y_fac;

close all;
figure()
ax1 = subplot(3,1,1);
plot(1:length(c1.y),c1.y);
hold on;
% for k3 = 1:k2
%     plot(OL_steps+k3:OL_steps+f+k3-1,c1.yfhat_k{k3})
% end
plot(OL_steps+1:OL_steps+sim_steps+f-1,r,'--')
xline(OL_steps+0.5,'k--');
ylabel('$y_k$','interpreter','latex');
grid on

ax2 = subplot(3,1,2);
plot(1:length(c1.u),c1.u);
hold on;
% for k3 = 1:k2
%     plot(OL_steps+k3:OL_steps+f+k3-1,c1.uf_k{k3})
% end
ylim([-16,16]);
xline(OL_steps+0.5,'k--');
yline(-15,'r--'); yline(15,'r--');
ylabel('$u_k$','interpreter','latex')
grid on

ax3 = subplot(3,1,3);
plot(1:length(e),e)
xline(OL_steps+0.5,'k--');
ylabel('$e_k$','interpreter','latex')
xlabel('samples','interpreter','latex')
grid on

linkaxes([ax1 ax2 ax3],'x')
xlim([1,k2+OL_steps])

function data = step_plant(data,e,plant,k1)
    data.y(:,k1)   = plant.C*data.x(:,k1) + plant.D*[data.u(:,k1); e(:,k1)];
    data.x(:,k1+1) = plant.A*data.x(:,k1) + plant.B*[data.u(:,k1); e(:,k1)];
end