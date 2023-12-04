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

addpath("../../data/raw",'-begin');
addpath(genpath("../../bin"),'-begin');
addpath(genpath("../../src"),'-begin');
rmpath(genpath("../../src/d_Re"));
rmpath(genpath("../../src/d_Nbar"));
addpath(genpath(pwd),'-begin');% add directory of this file to path

% cluster settings
cluster_flag = 2;
switch cluster_flag
    case 1 % single node
        myCluster = parcluster('local');
        Pool = parpool(myCluster, 4);
    case 2 % multiple nodes
        nworker = 120;
        myCluster = parcluster('SlurmProfile1')
        myCluster.ResourceTemplate = strjoin({'--job-name=d_pf', '--partition=compute',...
            '--time=17:00:00 --account=research-3mE-dcsc --nodes=10 --ntasks-per-node=12',... 17:00:00, 10, 12(/13?)
            '--cpus-per-task=1 --mem-per-cpu=4GB --output=d_pf.%j_%t.out --error=d_pf.%j_%t.err',.../dev/null
            '--mail-user=r.t.o.dinkla@tudelft.nl --mail-type=ALL'})
        Pool = parpool(myCluster, nworker);
end

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
p_min = 62;
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
Rdu= Ru/4;  % CL input disturbance variance

% OL-sim initial state
x0 = zeros(nx,1);

% number of CL simulation steps
CL_sim_steps = 1800;
num_steps = Nbar + CL_sim_steps;  % total simulation length

% define reference
ref = 50*[-sign(sin(2*pi/2*0.01*(0:CL_sim_steps-1))) ones(1,f_all(end)-1)]+50;

%% Running Simulations
% get output file name
if cluster_flag == 1
    files = dir(fullfile(pwd, 'd_pf.*.out'));
    fileNumbers = cellfun(@(x) str2double(regexp(x, 'd_pf\.(\d+)\.out', 'tokens', 'once')), {files.name});
    [~, maxIndex] = max(fileNumbers);
    outfile = files(maxIndex).name;
end

% description of data set
descr = strcat('Varying_pf_',num2str(p_min),'-',num2str(p_max),'-',num2str(num_p),...
    '_Nbar_',num2str(Nbar),'_Re_',num2str(Re),'_Ru_',num2str(Ru),'_Rdu_',num2str(Rdu),...
    '_Q_',num2str(Qk),'_R_',num2str(Rk),'_dR_',num2str(dRk));

% make directory for raw data files
raw_dir   = fullfile('..','..','data','raw' ,descr);
mkdir(raw_dir);

for k_p = 1:num_p
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
    mkdir(run_dir);
    tParFor = tic;
    % loop over noise realizations
    parfor k_e = 1:num_e
        seed_num = (k_p-1)*num_e+k_e+2520;
        loop_var(x0,N_OL,N_CL,p,f,k_p,k_e,plant,Ru,Re,ny,nu,nx,num_steps,Nbar,ref,Qk,Rk,dRk,num_c,Rdu,CL_sim_steps,run_dir,seed_num,Obsv_pf,Lu_pf,Ly_pf,Gu_pf);
    end
    dtParFor=toc(tParFor);
    clear output file contents
    if cluster_flag == 1
        fid = fopen(outfile,'w');
        fclose(fid);
    end

    % saved temp data into results structure and save in raw data folder
    temp2raw(num_e,num_c,run_dir,raw_dir,desc_runs,p,f,ref);
    disp(['End of run with p=f=',num2str(p),' Time taken=',num2str(dtParFor)])
    whos
end
disp('Done with iterations over p & f');
pause(1)
disp('Deleting parallel pool');
delete(Pool);
disp('Done');