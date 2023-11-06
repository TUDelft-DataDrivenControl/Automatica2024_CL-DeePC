rng default
clear all;
close all;
% clc;

%% Simulation settings
model_Favoreel1999 % loads model from Favoreel 1999 - original SPC paper

% controller settings
p = 20;
f = 20;
Ns_CL = (p+1)*nu + p*ny; %N >= n + (p+1)*nu + p*ny
Ns_OL = (p+f)*nu + p*ny; %N >= n + (p+f)*nu + p*ny
Nsbar_CL = p+Ns_CL;
Nsbar_OL = p+f+Ns_OL-1;
N_OL    = 500; %Nmin;
N_CL    = N_OL+f-1; % such that Nbar_CL = Nbar_OL
Nbar_CL = p+N_CL;
Nbar_OL = p+f+N_OL-1;
Qk = 100;
Rk = 0;
dRk= 10;

% number of controllers
num_c = 1;

% simulation length
OL_sim_steps = 1200;
CL_sim_steps = 500;
num_steps = OL_sim_steps + CL_sim_steps;

r = nan(ny,CL_sim_steps+f-1); % +f-1 needed for simulation end
r = (-square((0:size(r,2)-1)*2*pi/(200)))*50+1*50;

% noise
Re = 0.5*eye(ny);                       % variance
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

% normalize data
u_ol = u_ol./u_fac;
y_ol = y_ol./y_fac;
r = r./y_fac;

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
Rdu = 0.0;
du = mvnrnd(zeros(nu,1),Rdu,CL_sim_steps).';
du_max = 3.75;
du(abs(du)>du_max) = sign(du(abs(du)>du_max))*du_max; % prevent infeasibility
du = du./u_fac;

% create user-defined constraints
con = struct('Opti',cell(1,2));
for k = 1
con(k).Opti = casadi.Opti('conic');
con(k).uf   = con(k).Opti.variable(nu,f);
con(k).u0   = con(k).Opti.parameter(nu,1);
con(k).expr = {con(k).uf <=  15./u_fac,...
            con(k).uf >= -15./u_fac,...
            con(k).u0 - con(k).uf(:,1) <=  du_max./u_fac,...
            con(k).u0 - con(k).uf(:,1) >= -du_max./u_fac,...
            con(k).uf(:,1:end-1)-con(k).uf(:,2:end) <=  du_max./u_fac,...
            con(k).uf(:,1:end-1)-con(k).uf(:,2:end) >= -du_max./u_fac};
end

% initialize controllers
% 1) CL-DeePC, with IV
u1 = u_ol(:,end-Nbar_CL+1:end);
y1 = y_ol(:,end-Nbar_CL+1:end);
c1(1).controller = CL_DeePC(u1,y1,p,f,N_CL,Qk,Rk,dRk,constr=con(1));%ExplicitPredictor=true,use_IV=true,UseOptimizer=true,
c1(1).label = 'CL-DeePC, IV';
c1(1).color = '#005AB5';%'#1D3E23';

k1 = OL_sim_steps + 1;
k2 = 1;

kc = 1;
% first computed input
[c1(kc).uf_k{k2},c1(kc).yfhat_k{k2}] = c1(kc).controller.solve(rf=r(:,k2:k2+f-1));