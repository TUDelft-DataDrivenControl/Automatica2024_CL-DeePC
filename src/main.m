%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Closed-loop Data-Enabled Predictive Control
%           Authors: R. Dinkla, S.P. Mulders, T. Oomen, J.W. van Wingerden
%           Submission to Automatica 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
yalmip('clear');
clear;
% rng('default');
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

rho_max = max(abs(eig(A-K*C)),[],'all');
cond1_fac = -1/(2*log(rho_max)); % cond1_fac*log(N) < p

Nmax = 10^4;                      % maximum number of columns
pmin = ceil(cond1_fac*log(Nmax)); % such that above always true

% let minimum N be determined by regular DeePC
Nmin = pmin*(2*nu+ny); % p=f
num_N = 25;
N_all = Nmax; %round(logspace(log10(Nmin),log10(Nmax),num_N));
p = 25; % > pmin with Nmax
f = p;
Nol_bar_all = N_all+p+f-1;

% % controller settings
% p = 20;
% f = 20;
% Ns_CL = (p+1)*nu + p*ny; %N >= n + (p+1)*nu + p*ny
% Ns_OL = (p+f)*nu + p*ny; %N >= n + (p+f)*nu + p*ny
% Nsbar_CL = p+Ns_CL;
% Nsbar_OL = p+f+Ns_OL-1;
% N_OL    = 500; %Nmin;
% N_CL    = N_OL+f-1; % such that Nbar_CL = Nbar_OL
% Nbar_CL = p+N_CL;
% Nbar_OL = p+f+N_OL-1;
Qk = 100;
Rk = 0;
dRk= 10;

sims_per_N = 20;

% number of controllers
num_c = 2;

% initialize data structure
c1 = struct('u', cell(num_c,num_N,sims_per_N),);
fields = {'y','x','Cont','label'};
for i_field = 1:length(fields)
    c1(:).(fields{i_field}) = c1(:).u;
end

% noise
Re = 0.1*eye(ny);                        % variance

% for k_N = 1:num_N
% for k_r = 1:sims_per_N

% simulation length
OL_sim_steps = Nmax;
CL_sim_steps = 1800;
num_steps = OL_sim_steps + CL_sim_steps;

% predefine inputs, outputs, states, reference
for kc = 1:num_c
    c1(kc).u = nan(nu,num_steps);
    c1(kc).y = nan(ny,num_steps);
    c1(kc).x = nan(nx,num_steps+1);
end
r = nan(ny,CL_sim_steps+f-1); % +f-1 needed for simulation end
r = (-square((0:size(r,2)-1)*2*pi/(500)))*50+1*50;


e  = mvnrnd(zeros(ny,1),Re,num_steps).'; % trajectory realization

%% initial open loop simulation
x0 = zeros(nx,1); % initial state
Ru = 1*eye(nu);   % noise variance
u_ol = mvnrnd(zeros(nu,1),Ru,OL_sim_steps).';

% simulate system
[y_ol,~,x_ol] = lsim(plant,[u_ol;e(:,1:OL_sim_steps)],[],x0);
y_ol = y_ol.'; x_ol = x_ol.';
x_k = plant.A*x_ol(:,end) + plant.B*[u_ol(:,end); e(:,OL_sim_steps)];

% normalization factors
u_fac = max(abs(u_ol),[],2);
y_fac = max(abs(y_ol),[],2);

figure()
subplot(2,1,2)
plot(u_ol)
subplot(2,1,1)
plot(y_ol)

% normalize data
u_ol = u_ol./u_fac;
y_ol = y_ol./y_fac;
r = r./y_fac;

for kc = 1:num_c
    c1(kc).u(:,1:OL_sim_steps) = u_ol;
    c1(kc).y(:,1:OL_sim_steps) = y_ol;
    c1(kc).x(:,1:OL_sim_steps+1) = [x_ol x_k];
end

% plant normalization
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
% initialize data structures for controllers
for kc = 1:num_c
    c1(kc).uf_k    = cell(CL_sim_steps,1);
    c1(kc).yfhat_k = c1(kc).uf_k;
end

% initialize controllers
% 1) DeePC, with IV
u1 = u_ol(:,end-Nbar_OL+1:end);
y1 = y_ol(:,end-Nbar_OL+1:end);
c1(1).controller = DeePC(u1,y1,p,f,N_OL,Qk,Rk,dRk,use_IV=true);
c1(1).label = 'DeePC, IV';
c1(1).color = '#DC3220';

% 2) CL-DeePC, with IV
u2 = u_ol(:,end-Nbar_CL+1:end);
y2 = y_ol(:,end-Nbar_CL+1:end);
c1(2).controller = CL_DeePC(u2,y2,p,f,N_CL,Qk,Rk,dRk,use_IV=true);
c1(2).label = 'CL-DeePC, IV, explicit';
c1(2).color = '#005AB5';

for kc = 1:num_c
    tic
    % set counters
    k1 = OL_sim_steps + 1;
    k2 = 1;
    
    % first computed input
    [c1(kc).uf_k{k2},c1(kc).yfhat_k{k2}] = c1(kc).controller.solve(rf=r(:,k2:k2+f-1));
    c1(kc).u(:,k1) = c1(kc).uf_k{k2}(:,1) + du(:,k2);
    c1(kc) = step_plant(c1(kc),e,plant,k1);
    LHSt = c1(kc).controller.LHS(1:(p+f)*nu+c1(kc).controller.fid*nu,:);
    c1(kc).sv(:,k2) = svd(LHSt);
    c1(kc).rcond(:,k2) = rcond(LHSt);
    
    for k1 = OL_sim_steps+2:num_steps
        k2 = k2 + 1;
        
        [c1(kc).uf_k{k2},c1(kc).yfhat_k{k2}] = c1(kc).controller.step(c1(kc).u(:,k1-1),c1(kc).y(:,k1-1), rf=r(:,k2:k2+f-1)); 
        c1(kc).u(:,k1) = c1(kc).uf_k{k2}(:,1) + du(:,k2);
        c1(kc) = step_plant(c1(kc),e,plant,k1);
        LHSt = c1(kc).controller.LHS(1:(p+f)*nu+c1(kc).controller.fid*nu,:);
        c1(kc).sv(:,k2) = svd(LHSt);
        c1(kc).rcond(:,k2) = rcond(LHSt);
        
        % plot progression
        if rem(k2,100)==0
            clc;
            round(k2/CL_sim_steps*100,2)
            plot_all(c1,1,CL_sim_steps,OL_sim_steps,f,r,e,k2,u_fac,y_fac)
        end
    end
    toc
    er = c1(kc).y(:,OL_sim_steps+1:end)-r(:,1:CL_sim_steps);
    u  = c1(kc).u(:,OL_sim_steps+1:end);
    du = u-c1(kc).u(:,OL_sim_steps:end-1);
    c1(kc).cost = er(:).'*kron(speye(CL_sim_steps),Qk)*er(:)...
                + u(:).'*kron(speye(CL_sim_steps),Rk)*u(:)...
                + du(:).'*kron(speye(CL_sim_steps),dRk)*du(:);
end
%%
plot_all(c1,1,CL_sim_steps,OL_sim_steps,f,r,e,k2,u_fac,y_fac)

%%
function plot_all(c1,fignum,sim_steps,OL_steps,f,r,e,k2,u_fac,y_fac)
    
    figure(fignum);
    clf;

    % subplot with outputs & reference
    ax1 = subplot(2,1,1);
    xline(ax1,OL_steps+0.5,'k--','HandleVisibility','off');hold on;
    xline(ax1,OL_steps+c1(1).controller.Nbar+0.5,'k--','HandleVisibility','off'); hold on;
    r = r.*y_fac; % reverse normalization of reference
    plot(ax1,OL_steps+1:OL_steps+sim_steps+f-1,r,'--',...
        'DisplayName','reference',...
        'color', [.5 .5 .5], ... grey
        'linewidth', 1.5);%      thicker line
    ylabel(ax1,'$y_k$','interpreter','latex');
    grid on
    
    % subplot with inputs
    ax2 = subplot(2,1,2);
    xline(ax2,OL_steps+0.5,'k--');hold on;
    xline(ax2,OL_steps+c1(1).controller.Nbar+0.5,'k--');
    yline(ax2,-15,'r--');
    yline(ax2, 15,'r--');
    ylim(ax2,[-16,16]);
    ylabel(ax2,'$u_k$','interpreter','latex')
    grid on
    
    num_c = numel(c1);
    lineref1 = cell(num_c,1);
    lineref2 = cell(num_c,1);
    Styles = struct('LineStyle',cell(num_c,1),'Color',cell(num_c,1),'LineWidth',cell(num_c,1));
    legend_descr = cell(num_c,1);
    for kc = num_c:-1:1
        legend_descr{kc} = c1(kc).label;
        % reverse normalization
        c1(kc).u = c1(kc).u.*u_fac;
        c1(kc).y = c1(kc).y.*y_fac;
        
        lineref1{kc} = plot(ax1,1:length(c1(kc).y),c1(kc).y,'DisplayName', c1(kc).label);
        lineref1{kc}.LineWidth = 1;
        lineref1{kc}.Color = c1(kc).color;
        Styles(kc).LineStyle = lineref1{kc}.LineStyle;
        Styles(kc).Color     = lineref1{kc}.Color;
        Styles(kc).LineWidth = lineref1{kc}.LineWidth;
        % for k3 = 1:k2
        %     plot(ax1,OL_steps+k3:OL_steps+f+k3-1,cont_struct.yfhat_k{k3})
        % end
        styles2use = namedargs2cell(Styles(kc));
        lineref2{kc} = plot(ax2,1:length(c1(kc).u),c1(kc).u,styles2use{:});
        % for k3 = 1:k2
        %     plot(ax2,OL_steps+k3:OL_steps+f+k3-1,cont_struct.uf_k{k3})
        % end
    end
    legend(ax1);

%     ax3 = subplot(3,1,3);
%     plot(ax3,1:length(e),e)
%     xline(ax3,OL_steps+0.5,'k--');
%     ylabel('$e_k$','interpreter','latex')
%     xlabel('samples','interpreter','latex')
%     grid on
    
    linkaxes([ax1 ax2],'x')
    xlim([OL_steps,k2+OL_steps])
end

function data = step_plant(data,e,plant,k1)
    data.y(:,k1)   = plant.C*data.x(:,k1) + plant.D*[data.u(:,k1); e(:,k1)];
    data.x(:,k1+1) = plant.A*data.x(:,k1) + plant.B*[data.u(:,k1); e(:,k1)];
end