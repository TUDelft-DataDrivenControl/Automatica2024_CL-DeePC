%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Closed-loop Data-Enabled Predictive Control
%           Authors: R. Dinkla, S.P. Mulders, T. Oomen, J.W. van Wingerden
%           Submission to Automatica 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear;
rng('default');
clc;

%% Edit path - add dependencies
[mfilePath,~,~] = fileparts(mfilename('fullpath')); % get dir of this file
cd(mfilePath); % change cwd to folder containing this file

addpath(genpath("../data"),'-begin')
addpath(genpath("../bin"),'-begin')
addpath(genpath(pwd),'-begin');% add directory of this file to path

%% Simulation settings
model_Favoreel1999 % loads model from Favoreel 1999 - original SPC paper

rho_max = max(abs(eig(A-K*C)),[],'all');
cond1_fac = -1/(2*log(rho_max)); % cond1_fac*log(N) < p

Nmax = 10^3;                      % maximum number of columns
pmin = ceil(cond1_fac*log(Nmax)); % such that above always true

% controller settings
p = 20; % > pmin with Nmax
f = p;
Nmin = p*(nu+ny)+f*nu; % minimum N is determined by regular DeePC
Qk = 100;
Rk = 0;
dRk= 10;

% number of
num_c = 2;   % controllers
num_e = 100; % noise realizations per value of N
num_N = 50;  % values for N

% N & Nbar values - same Nbar for DeePC & CL-DeePC
N_all_OL = round(logspace(log10(Nmin),log10(Nmax),num_N));
N_all_CL = N_all_OL + f-1;
Nbar_all = N_all_OL+p+f-1;
Nbar_min = min(Nbar_all,[],'all');
Nbar_max = max(Nbar_all,[],'all');

results.Nbar  = Nbar_all;
results.CzLabel = {'DeePC, IV','CL-DeePC, IV'};           % labels
results.noise = cell(num_N,num_e);       % innovation noise
results.du_CL = cell(num_N,num_e);
results.u_OL  = cell(num_N,num_e);
results.y_OL  = cell(num_N,num_e);
results.x_OL  = cell(num_N,num_e);
results.u_CL  = cell(num_N,num_e,num_c); % CL inputs
results.y_CL  = cell(num_N,num_e,num_c); % CL outputs
results.x_CL  = cell(num_N,num_e,num_c); % CL states
results.Cost  = cell(num_N,num_e,num_c);

% variances
Re = 0.25*eye(ny); % noise
Ru = 1*eye(nu);   % OL input
Rdu= 0.01;       % CL input disturbance

% OL-sim initial state
x0 = zeros(nx,1);

% number of CL simulation steps
CL_sim_steps = 1800;

% define reference
ref = nan(ny,CL_sim_steps+f-1); % +f-1 needed for simulation end
ref = (-square((0:size(ref,2)-1)*2*pi/(200)))*50+1*50;

%% Run Simulations
temp_str = 'Varying_Nbar_kn_';
parpool(10);
for k_N = 1:num_N
Nbar = Nbar_all(k_N);
N_OL = N_all_OL(k_N);
N_CL = N_all_CL(k_N);
num_steps = Nbar + CL_sim_steps;  % total simulation length

    parfor k_e = 1:num_e
        loop_var(x0,N_OL,N_CL,p,f,k_N,k_e,plant,Ru,Re,ny,nu,nx,num_steps,Nbar,ref,Qk,Rk,dRk,num_c,Rdu,CL_sim_steps,temp_str);
    end

        % put saved temporary data into results structure
    for k_e = 1:num_e
        Str = strcat('../data/temp/',temp_str,num2str(k_N),'_ke_',num2str(k_e),'.mat');
        clear u_OL y_OL x_OL e du_CL u_CL y_CL x_CL Cost
        load(Str,'u_OL','y_OL','x_OL','e','du_CL','u_CL','y_CL','x_CL','Cost');
        results.noise{k_N,k_e} = e;
        results.du_CL{k_N,k_e} = du_CL;
        results.u_OL{k_N,k_e}  = u_OL;
        results.y_OL{k_N,k_e}  = y_OL;
        results.x_OL{k_N,k_e}  = x_OL;
        results.u_CL(k_N,k_e,:)  = u_CL; % CL inputs
        results.y_CL(k_N,k_e,:)  = y_CL; % CL outputs
        results.x_CL(k_N,k_e,:)  = x_CL; % CL states
        results.Cost(k_N,k_e,:)  = Cost;
        delete(Str);
    end

    % save results structure
    save(strcat('../data/raw/Varying_Nbar/Varying_Nbar_',num2str(Nbar_min),'-',num2str(Nbar_max),'-',num2str(num_N),...
        '_p_',num2str(p),'_f_',num2str(f),'_eVar_',num2str(Re),'_uPEVar_',num2str(Rdu),'_Q_',num2str(Qk),'_R_',...
        num2str(Rk),'_dR_',num2str(dRk),'.mat'),'results');

end