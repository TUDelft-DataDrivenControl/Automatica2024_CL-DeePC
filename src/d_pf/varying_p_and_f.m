function varying_p_and_f(k_p,k_e)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Closed-loop Data-Enabled Predictive Control
%           Authors: R. Dinkla, S.P. Mulders, T. Oomen, J.W. van Wingerden
%           Submission to Automatica 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
%clear;
rng('default');
clc;

%% Edit path - add dependencies
[mfilePath,~,~] = fileparts(mfilename('fullpath')); % get dir of this file
cd(mfilePath); % change cwd to folder containing this file

addpath(fullfile("..","..","data","raw"),'-begin');
addpath(genpath(fullfile("..","..","bin")),'-begin');
addpath(genpath(fullfile("..","..","src")),'-begin');
rmpath(genpath(fullfile("..","..","src","d_Re")));
rmpath(genpath(fullfile("..","..","src","d_Nbar")));
rmpath(genpath(fullfile("..","..","src","test")));
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
num_e = 120; % noise realizations per value of p,f

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

[Obsv_f,Lu_act,Ly_act,Gu_act] = make_ObsvLuLyGu(plant,f_all(end),p_all(end));

% variances
Re = 0.25*eye(ny); % noise
Ru = 1*eye(nu);    % OL input
Rdu= 0*Ru/4;  % CL input disturbance variance

% OL-sim initial state
x0 = zeros(nx,1);

% number of CL simulation steps
CL_sim_steps = 1800;
num_steps = Nbar + CL_sim_steps;  % total simulation length

% define reference
ref = 50*[-sign(sin(2*pi/2*0.01*(0:CL_sim_steps-1))) ones(1,f_all(end)-1)]+50;

%% Running Simulation

% description of data set
descr = strcat('Varying_pf_',num2str(p_min),'-',num2str(p_max),'-',num2str(num_p),...
    '_Nbar_',num2str(Nbar),'_Re_',num2str(Re),'_Ru_',num2str(Ru),'_Rdu_',num2str(Rdu),...
    '_Q_',num2str(Qk),'_R_',num2str(Rk),'_dR_',num2str(dRk));

% make directory for raw data files
raw_dir   = fullfile('..','..','data','raw' ,descr);
if ~isfolder(raw_dir)
    mkdir(raw_dir);
end

N_OL = N_OL_all(k_p);
N_CL = N_CL_all(k_p);
p = p_all(k_p);
f = f_all(k_p);

Obsv_pf = Obsv_f(1:ny*f,:);
Lu_pf   = Lu_act(1:ny*f,1:nu*p);
Ly_pf   = Ly_act(1:ny*f,1:ny*p);
Gu_pf   = Gu_act(1:ny*f,1:nu*f);

% make directories to save data in
desc_runs = strcat('p_',num2str(p),'_f_',num2str(f));
run_dir   = fullfile(raw_dir,desc_runs);
if ~isfolder(run_dir)
    mkdir(run_dir); addpath(run_dir);
end

tic
seed_num = (k_p-1)*num_e+k_e;
loop_var(x0,N_OL,N_CL,p,f,k_p,k_e,plant,Ru,Re,ny,nu,nx,num_steps,Nbar,ref,Qk,Rk,dRk,num_c,Rdu,CL_sim_steps,run_dir,seed_num,Obsv_pf,Lu_pf,Ly_pf,Gu_pf);
toc

% saved temp data into results structure and save in raw data folder
%temp2raw(num_e,num_c,run_dir,raw_dir,desc_runs,p,f,ref);
end