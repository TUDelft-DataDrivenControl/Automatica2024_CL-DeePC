%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Closed-loop Data-Enabled Predictive Control
%           Authors: R. Dinkla, S.P. Mulders, T. Oomen, J.W. van Wingerden
%           Submission to Automatica 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear;
% rng default;
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
for k=1:length(dirs)
    addpath(genpath(dirs{k}),'-begin')
end

% add directory of this file to path
addpath(genpath(pwd),'-begin');

%% Simulation settings
model_Favoreel1999 % loads model from Favoreel 1999 - original SPC paper

% controller settings
p = 15;
f = 50;
Nmin = (p+1)*nu + p*ny; %N >= n + (p+1)*nu + p*ny
N  = 1*Nmin;
Nbar = p+N;
Qk = 10;
Rk = 0.01;

% simulation length
num_steps = Nbar*5;

% predefine inputs, outputs, states, reference
u = nan(nu,num_steps);
y = nan(ny,num_steps);
x = nan(nx,num_steps);
r = nan(ny,num_steps+f-1); % +f-1 for simulation end
r = (-square((0:size(r,2)-1)*2*pi/Nbar))*50-0.5*50;

% noise
Re = 0.1*eye(ny);                        % variance
e  = mvnrnd(zeros(ny,1),Re,num_steps).'; % trajectory realization

%% initial open loop simulation
x0 = zeros(nx,1); % initial state
Ru = 1*eye(nu);   % noise variance
u(:,1:Nbar)  = mvnrnd(zeros(nu,1),Ru,Nbar).';

% simulate system
[y_ol,~,x_ol] = lsim(plant,[u(:,1:Nbar);e(:,1:Nbar)],[],x0);
y(:,1:Nbar) = y_ol.';
x(:,1:Nbar) = x_ol.';

%% closed-loop operation

% initialize controllers
Cont1 = CL_DeePC(u(:,1:Nbar),y(:,1:Nbar),p,f,N,Qk,Rk);

[uf,yf_hat,G] = Cont1.optimize_solve(r(:,Nbar+1:Nbar+f));

%%
close all;
%%
figure()
ax1 = subplot(2,1,1);
plot(u);
hold on;
plot(Nbar+1:Nbar+f,uf)
grid on

ax2 = subplot(2,1,2);
plot(y);
hold on;
plot(Nbar+1:Nbar+f,yf_hat)
plot(r(:,1:Nbar+f),'--')
grid on

linkaxes([ax1 ax2],'x')