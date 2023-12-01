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

addpath(genpath("../../data"),'-begin')
addpath(genpath("../../bin"),'-begin')
addpath(genpath("../../src"),'-begin')
addpath(genpath(pwd),'-begin');% add directory of this file to path

% cluster settings
% myCluster = parcluster('local');
% parpool(myCluster, 4);
%nworker = 48;
%myCluster = parcluster('SlurmProfile1')
%myCluster.ResourceTemplate = strjoin({'--job-name=d_Nbar', '--partition=compute',...
%    '--time=08:00:00 --account=research-3mE-dcsc --nodes=4 --ntasks=48',...
%    '--cpus-per-task=1 --mem-per-cpu=4GB --output=d_Nbar.%j.out --error=d_Nbar.%j.err',...
%    '--mail-user=r.t.o.dinkla@tudelft.nl --mail-type=ALL'})
%parpool(myCluster, nworker)

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

[Obsv_f,Lu_act,Ly_act,Gu_act] = make_ObsvLuLyGu(plant,f,p);

% number of
num_c = 2;   % controllers
num_e = 4;%120; % noise realizations per value of N
num_N = 1;%50;  % values for N

% N & Nbar values - same Nbar for DeePC & CL-DeePC
N_all_OL = Nmax;%round(logspace(log10(Nmin),log10(Nmax),num_N));
N_all_CL = N_all_OL + f-1;
Nbar_all = N_all_OL+p+f-1;
Nbar_min = min(Nbar_all,[],'all');
Nbar_max = max(Nbar_all,[],'all');

% variances
Re = 0.25*eye(ny); % noise
Ru = 1*eye(nu);    % OL input
Rdu= Ru/4;         % CL input disturbance

% OL-sim initial state
x0 = zeros(nx,1);

% number of CL simulation steps
OL_sim_steps = 1200; %>=Nbar_max
CL_sim_steps = 1800;
num_steps = OL_sim_steps + CL_sim_steps;  % total simulation length

% define reference
ref = 50*[-sign(sin(2*pi/2*0.01*(0:CL_sim_steps-1))) ones(1,f-1)]+50;

%% Run Simulations
% get output file name
% files = dir(fullfile(pwd, 'd_Nbar.*.out'));
% fileNumbers = cellfun(@(x) str2double(regexp(x, 'd_Nbar\.(\d+)\.out', 'tokens', 'once')), {files.name});
% [~, maxIndex] = max(fileNumbers);
% outfile = files(maxIndex).name;

% description of data set
descr = strcat('Varying_Nbar_',num2str(Nbar_min),'-',num2str(Nbar_max),'-',num2str(num_N),...
    '_p_',num2str(p), '_f_',num2str(f), '_Re_',num2str(Re),'_Ru_',num2str(Ru),'_Rdu_',num2str(Rdu),...
    '_Q_',num2str(Qk),'_R_',num2str(Rk),'_dR_',num2str(dRk));

% make directory for raw data files
raw_dir   = fullfile('..','..','data','raw' ,descr);
mkdir(raw_dir);

tic
for k_N = 1:num_N
    Nbar = Nbar_all(k_N);
    N_OL = N_all_OL(k_N);
    N_CL = N_all_CL(k_N);
    
    % make directories to save data in
    desc_runs = strcat('Nbar_',num2str(Nbar));
    run_dir   = fullfile(raw_dir,desc_runs);
    mkdir(run_dir);
    
    % loop over noise realizations
    parfor k_e = 1:num_e
        seed_num = (k_N-1)*num_e+k_e;
        loop_var(x0,N_OL,N_CL,p,f,k_N,k_e,plant,Ru,Re,ny,nu,nx,num_steps,Nbar,ref,Qk,Rk,dRk,num_c,Rdu,CL_sim_steps,run_dir,seed_num,Obsv_f,Lu_act,Ly_act,Gu_act);
    end

    % clear output file contents
%     fid = fopen(outfile,'w');
%     fclose(fid);

    % saved temp data into results structure and save in raw data folder
    temp2raw(num_e,num_c,run_dir,raw_dir,desc_runs,Nbar,ref);
end
toc