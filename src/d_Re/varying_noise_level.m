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
addpath(genpath(pwd),'-begin');% add directory of this file to path

% cluster settings
myCluster = parcluster('local');
parpool(myCluster, 41);
%nworker = 48;
%myCluster = parcluster('SlurmProfile1')
%myCluster.ResourceTemplate = strjoin({'--job-name=d_Re', '--partition=compute',...
%    '--time=08:00:00 --account=research-3mE-dcsc --nodes=4 --ntasks=48',...
%    '--cpus-per-task=1 --mem-per-cpu=4GB --output=d_Re.%j.out --error=d_Re.%j.err',...
%    '--mail-user=r.t.o.dinkla@tudelft.nl --mail-type=ALL'})
%parpool(myCluster, nworker)

%% Simulation settings
model_Favoreel1999 % loads model from Favoreel 1999 - original SPC paper

rho_max = max(abs(eig(A-K*C)),[],'all');
cond1_fac = -1/(2*log(rho_max)); % cond1_fac*log(N) < p

% controller settings
p = 20;
f = p;
Nmin = p*(nu+ny)+f*nu; % minimum N is determined by regular DeePC
N_OL = 200;
N_CL = N_OL + f-1;     % s.t. Nbar remains the same
Nbar = N_OL + p + f-1; % = N_CL + p;
pmin = ceil(cond1_fac*log(N_CL)); % check that p > pmin
Qk = 100;
Rk = 0;
dRk= 10;

[Obsv_f,Lu_act,Ly_act,Gu_act] = make_ObsvLuLyGu(plant,f,p);

% number of
num_c = 2;     % controllers
num_n = 50;     % number of noise levels
num_e = 120; % number of noise realizations

% variances
Re_all = logspace(-4,0,num_n);
Re_min = min(Re_all,[],'all');
Re_max = max(Re_all,[],'all');
Ru = 1*eye(nu);   % OL input
Rdu= Ru/4;       % CL input disturbance

% OL-sim initial state
x0 = zeros(nx,1);

% number of simulation steps
CL_sim_steps = 1800;
num_steps = Nbar + CL_sim_steps;  % total simulation length

% define reference
ref = 50*[-sign(sin(2*pi/2*0.01*(0:CL_sim_steps-1))) ones(1,f-1)]+50;

%% running simulations
% get output file name
files = dir(fullfile(pwd, 'd_Re.*.out'));
fileNumbers = cellfun(@(x) str2double(regexp(x, 'd_Re\.(\d+)\.out', 'tokens', 'once')), {files.name});
[~, maxIndex] = max(fileNumbers);
outfile = files(maxIndex).name;

% description of data set
descr = strcat('Varying_Re_',num2str(Re_min),'-',num2str(Re_max),'-',num2str(num_n),...
    '_Nbar_',num2str(Nbar),'_p_',num2str(p),'_f_',num2str(f),'_Ru_',num2str(Ru),'_Rdu_',num2str(Rdu),...
    '_Q_',num2str(Qk),'_R_',num2str(Rk),'_dR_',num2str(dRk));

% make directory for raw data files
raw_dir   = fullfile('..','..','data','raw' ,descr);
mkdir(raw_dir);

for k_n = 1:num_n % loop over noise levels
    Re = Re_all(k_n)*eye(ny);

    % make directories to save data in
    desc_runs = strcat('Re_',num2str(Re));
    run_dir   = fullfile(raw_dir,desc_runs);
    mkdir(run_dir);

    % loop over noise realizations
    parfor k_e = 1:num_e
        seed_num = (k_n-1)*num_e+k_e;
        loop_var(x0,N_OL,N_CL,p,f,k_n,k_e,plant,Ru,Re,ny,nu,nx,num_steps,Nbar,ref,Qk,Rk,dRk,num_c,Rdu,CL_sim_steps,run_dir,seed_num,Obsv_f,Lu_act,Ly_act,Gu_act);
    end

    % clear output file contents
    fid = fopen(outfile,'w');
    fclose(fid);
end