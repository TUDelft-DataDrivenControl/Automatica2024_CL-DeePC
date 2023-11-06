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
Plant = plant;

rho_max = max(abs(eig(A-K*C)),[],'all');
cond1_fac = -1/(2*log(rho_max)); % cond1_fac*log(N) < p

% controller settings
p = 20;
f = p;
Nmin = p*(nu+ny)+f*nu; % minimum N is determined by regular DeePC
N_OL = 500;
N_CL = N_OL + f-1;     % s.t. Nbar remains the same
Nbar = N_OL + p + f-1; % = N_CL + p;
pmin = ceil(cond1_fac*log(N_CL)); % check that p > pmin
Qk = 100;
Rk = 0;
dRk= 10;

% number of
num_c = 2;     % controllers
num_n = 50;     % number of noise levels
num_e = 100; % number of noise realizations

% variances
eVar = logspace(-4,0,num_n);
Re_min = min(eVar,[],'all');
Re_max = max(eVar,[],'all');
Ru = 1*eye(nu);   % OL input
Rdu= 0.01;       % CL input disturbance

% N & Nbar values - same Nbar for DeePC & CL-DeePC

% initialize data cells
results.eVar  = eVar;
results.noise = cell(num_n,num_e);       % innovation noise
results.du_CL = cell(num_n,num_e);
results.u_OL  = cell(num_n,num_e);
results.y_OL  = cell(num_n,num_e);
results.x_OL  = cell(num_n,num_e);
results.u_CL  = cell(num_n,num_e,num_c); % CL inputs
results.y_CL  = cell(num_n,num_e,num_c); % CL outputs
results.x_CL  = cell(num_n,num_e,num_c); % CL states
results.Cost  = cell(num_n,num_e,num_c);
results.CzLabel = {'DeePC, IV','CL-DeePC, IV'};           % labels

% OL-sim initial state
x0 = zeros(nx,1);

% number of simulation steps
CL_sim_steps = 1800;
num_steps = Nbar + CL_sim_steps;  % total simulation length

% define reference
ref = nan(ny,CL_sim_steps+f-1); % +f-1 needed for simulation end
ref = (-square((0:size(ref,2)-1)*2*pi/(200)))*50+1*50;


%% running simulations
temp_str = 'Varying_eVar_kn_';
parpool(10);
for k_n = 1:num_n % loop over noise levels
    Re = eVar(k_n)*eye(ny);
    parfor k_e = 1:num_e % loop over noise realizations
        loop_var(x0,N_OL,N_CL,p,f,k_n,k_e,plant,Ru,Re,ny,nu,nx,num_steps,Nbar,ref,Qk,Rk,dRk,num_c,Rdu,CL_sim_steps,temp_str);
    end

    % put saved temporary data into results structure
    for k_e = 1:num_e
        Str = strcat('../data/temp/',temp_str,num2str(k_n),'_ke_',num2str(k_e),'.mat');
        clear u_OL y_OL x_OL e du_CL u_CL y_CL x_CL Cost
        load(Str,'u_OL','y_OL','x_OL','e','du_CL','u_CL','y_CL','x_CL','Cost');
        results.noise{k_n,k_e} = e;
        results.du_CL{k_n,k_e} = du_CL;
        results.u_OL{k_n,k_e}  = u_OL;
        results.y_OL{k_n,k_e}  = y_OL;
        results.x_OL{k_n,k_e}  = x_OL;
        results.u_CL(k_n,k_e,:)  = u_CL; % CL inputs
        results.y_CL(k_n,k_e,:)  = y_CL; % CL outputs
        results.x_CL(k_n,k_e,:)  = x_CL; % CL states
        results.Cost(k_n,k_e,:)  = Cost;
        delete(Str);
    end

    % save results structure
    save(strcat('../data/raw/Varying_eVar/Varying_eVar_',num2str(Re_min),'-',num2str(Re_max),'-',num2str(num_n),...
    '_Nbar_',num2str(Nbar),'_p_',num2str(p),'_f_',num2str(f),'_uPEVar_',num2str(Rdu),...
    '_Q_',num2str(Qk),'_R_',num2str(Rk),'_dR_',num2str(dRk),'.mat'),'results');
end