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

% N & Nbar values - same Nbar for DeePC & CL-DeePC

% initialize data cells
results.eVar  = Re_all;
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
ref = 50*[-sign(sin(2*pi/2*0.01*(0:CL_sim_steps-1))) ones(1,f-1)]+50;

%% running simulations
% get output file name
files = dir(fullfile(pwd, 'd_Re.*.out'));
fileNumbers = cellfun(@(x) str2double(regexp(x, 'd_Re\.(\d+)\.out', 'tokens', 'once')), {files.name});
[~, maxIndex] = max(fileNumbers);
outfile = files(maxIndex).name;

myCluster = parcluster('local');
parpool(myCluster, 41);
descr = strcat('Varying_Re_',num2str(Re_min),'-',num2str(Re_max),'-',num2str(num_n),...
    '_Nbar_',num2str(Nbar),'_p_',num2str(p),'_f_',num2str(f),'_Ru_',num2str(Ru),'_Rdu_',num2str(Rdu),...
    '_Q_',num2str(Qk),'_R_',num2str(Rk),'_dR_',num2str(dRk));
for k_n = 1:num_n % loop over noise levels
    Re = Re_all(k_n)*eye(ny);

    % make directories to save data in
    desc_runs = strcat('Re_',num2str(Re));
    temp_dir1 = fullfile('..','data','temp',descr);
    temp_dir2 = fullfile(temp_dir1,desc_runs);
    raw_dir   = fullfile('..','data','raw' ,descr);
    mkdir(temp_dir2);
    mkdir(raw_dir)

    % loop over noise realizations
    parfor k_e = 1:num_e
        seed_num = (num_n-1)*num_e+k_e;
        loop_var(x0,N_OL,N_CL,p,f,k_n,k_e,plant,Ru,Re,ny,nu,nx,num_steps,Nbar,ref,Qk,Rk,dRk,num_c,Rdu,CL_sim_steps,temp_dir2,seed_num,Obsv_f,Lu_act,Ly_act,Gu_act);
    end

    % clear output file contents
    fid = fopen(outfile,'w');
    fclose(fid);
    
    % put saved temporary data into results structure
    temp2raw(num_e,num_c,temp_dir1,temp_dir2,raw_dir,desc_runs,Re);
end