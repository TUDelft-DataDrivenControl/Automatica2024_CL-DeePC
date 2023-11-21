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

% controller settings
Qk = 100;
Rk = 0;
dRk= 10;

% number of
num_c = 2; % controllers
num_e = 100; % noise realizations per value of p,f

% p & f values
p_min = 20;
p_max = 100;
p_all = p_min:2:p_max; % > pmin with Nmax
num_p = length(p_all);  % number of values for p
f_all = p_all; % p = f
Nmin = p_all*(nu+ny)+p_all*nu;
Nmax = floor(rho_max.^(-2*p_all));

Nbar = 10^3 + 2*p_max -1; % for both DeePC & CL-DeePC
N_CL_all = Nbar - p_all;
N_OL_all = Nbar - p_all - f_all + 1;

% initialize data cells
results.p     = p_all;
results.f     = f_all;
results.CzLabel = {'DeePC, IV','CL-DeePC, IV'};           % labels
results.noise = cell(num_p,num_e);       % innovation noise
results.du_CL = cell(num_p,num_e);
results.u_OL  = cell(num_p,num_e);
results.y_OL  = cell(num_p,num_e);
results.x_OL  = cell(num_p,num_e);
results.u_CL  = cell(num_p,num_e,num_c); % CL inputs
results.y_CL  = cell(num_p,num_e,num_c); % CL outputs
results.x_CL  = cell(num_p,num_e,num_c); % CL states
results.Cost  = cell(num_p,num_e,num_c); % cost of controller implementation

% variances
Re = 0.25*eye(ny); % noise
Ru = 0.1*eye(nu);    % OL input
Rdu= Ru/4;  % CL input disturbance variance

% OL-sim initial state
x0 = zeros(nx,1);

% number of CL simulation steps
CL_sim_steps = 1800;
num_steps = Nbar + CL_sim_steps;  % total simulation length

% define reference
ref = 50*[-sign(sin(2*pi/2*0.01*(0:num_steps-1))) ones(1,f_all(end))]+50;

%% Running Simulations
temp_str = 'Varying_p_kp_';
parpool(10);
seed_nums = reshape(1:num_p*num_e,num_p,num_e);
for k_p = 1:num_p
    N_OL = N_OL_all(k_p);
    N_CL = N_CL_all(k_p);
    p = p_all(k_p);
    f = f_all(k_p);
    
    % loop over noise realizations
    parfor k_e = 1:num_e
        seed_num = (num_p-1)*num_e+k_e;
        loop_var(x0,N_OL,N_CL,p,f,k_p,k_e,plant,Ru,Re,ny,nu,nx,num_steps,Nbar,ref,Qk,Rk,dRk,num_c,Rdu,CL_sim_steps,temp_str,seed_num);
    end
    
    % put saved temporary data into results structure
    for k_e = 1:num_e
        Str = strcat('../data/temp/',temp_str,num2str(k_p),'_ke_',num2str(k_e),'.mat');
        clear u_OL y_OL x_OL e du_CL u_CL y_CL x_CL Cost
        load(Str,'u_OL','y_OL','x_OL','e','du_CL','u_CL','y_CL','x_CL','Cost');
        results.noise{k_p,k_e} = e;
        results.du_CL{k_p,k_e} = du_CL;
        results.u_OL{k_p,k_e}  = u_OL;
        results.y_OL{k_p,k_e}  = y_OL;
        results.x_OL{k_p,k_e}  = x_OL;
        results.u_CL(k_p,k_e,:)  = u_CL; % CL inputs
        results.y_CL(k_p,k_e,:)  = y_CL; % CL outputs
        results.x_CL(k_p,k_e,:)  = x_CL; % CL states
        results.Cost(k_p,k_e,:)  = Cost;
        delete(Str);
    end

    % save results structure
    save(strcat('../data/raw/Varying_p/Varying_p_',num2str(p_min),'-',num2str(p_max),'-',num2str(num_p),...
    '_Nbar_',num2str(Nbar),'_eVar_',num2str(Re),'_uPEVar_',num2str(Rdu),...
    '_Q_',num2str(Qk),'_R_',num2str(Rk),'_dR_',num2str(dRk),'.mat'),'results');
end