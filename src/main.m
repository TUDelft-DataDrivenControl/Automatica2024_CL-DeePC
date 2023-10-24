%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Closed-loop Data-Enabled Predictive Control
%           Authors: R. Dinkla, S.P. Mulders, T. Oomen, J.W. van Wingerden
%           Submission to Automatica 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
yalmip('clear');
clear;
% rng('default');
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
for k1=1:length(dirs)
    addpath(genpath(dirs{k1}),'-begin')
end

% add directory of this file to path
addpath(genpath(pwd),'-begin');

%% Simulation settings
model_Favoreel1999 % loads model from Favoreel 1999 - original SPC paper
Plant = plant;

rho_max = max(abs(eig(A-K*C)),[],'all');
cond1_fac = -1/(2*log(rho_max)); % cond1_fac*log(N) < p

Nmax = 10^3;                      % maximum number of columns
pmin = ceil(cond1_fac*log(Nmax)); % such that above always true

% minimum N is determined by regular DeePC
p = 25; % > pmin with Nmax
f = p;
Nmin = p*(nu+ny)+f*nu; % p=f

% % controller settings
% p = 20;
% f = 20;
% Ns_CL = (p+1)*nu + p*ny; %N >= n + (p+1)*nu + p*ny
% Ns_OL = (p+f)*nu + p*ny; %N >= n + (p+f)*nu + p*ny
% Nsbar_CL = p+Ns_CL;
% Nsbar_OL = p+f+Ns_OL-1;
% N_OL    = 500; %Nmin;
% N_CL    = N_OL+f-1; % such that Nbar_CL = Nbar_OL
% Nbar_CL = p+N_CL;
% Nbar_OL = p+f+N_OL-1;
Qk = 100;
Rk = 0;
dRk= 10;

% number of
num_c = 2;   % controllers
num_e = 100; % noise realizations per value of N
num_N = 25;  % values for N

% N & Nbar values - same Nbar for DeePC & CL-DeePC
N_all_OL = round(logspace(log10(Nmin),log10(Nmax),num_N));
N_all_CL = N_all_OL + f-1;
Nbar_all = N_all_OL+p+f-1;

% initialize data structure
% fields = {'u','y','x','e','r','Cont','label'};
% data = cell(num_N, num_c, num_e, length(fields));
% data = cell2struct(data,fields,4);

noise = cell(num_N,num_e);       % innovation noise
du_CL = cell(num_N,num_e);
u_OL  = cell(num_N,num_e);
y_OL  = cell(num_N,num_e);
x_OL  = cell(num_N,num_e);
Cz    = cell(num_N,num_e,num_c); % controllers
Label = {'DeePC, IV','CL-DeePC, IV'};           % labels
Color = {'#DC3220','#005AB5'};
u_CL  = cell(num_N,num_e,num_c); % CL inputs
y_CL  = cell(num_N,num_e,num_c); % CL outputs
x_CL  = cell(num_N,num_e,num_c); % CL states
Cost  = cell(num_N,num_e,num_c);

% variances
Re = 0.1*eye(ny); % noise
Ru = 1*eye(nu);   % OL input
Rdu= Ru/10;       % CL input disturbance

% OL-sim initial state
x0 = zeros(nx,1);

% number of CL simulation steps
CL_sim_steps = 2000;

% define reference
ref = nan(ny,CL_sim_steps+f-1); % +f-1 needed for simulation end
ref = (-square((0:size(ref,2)-1)*2*pi/(500)))*50+1*50;

for k_N = 1%1:num_N
    % total simulation length
Nbar = Nbar_all(k_N);
N_OL = N_all_OL(k_N);
N_CL = N_all_CL(k_N);
num_steps = Nbar + CL_sim_steps;

tic
for k_e = 1%:num_e
plant = Plant;

% noise realization
e  = mvnrnd(zeros(ny,1),Re,num_steps).';
noise{k_N,k_e} = e;

%% initial open loop simulation
u_ol = mvnrnd(zeros(nu,1),Ru,Nbar).';

% simulate system
[y_ol,~,x_ol] = lsim(plant,[u_ol;e(:,1:Nbar)],[],x0);
y_ol = y_ol.'; x_ol = x_ol.';
x0_CL = plant.A*x_ol(:,end) + plant.B*[u_ol(:,end); e(:,Nbar)];

% saving data
u_OL{k_N,k_e} = u_ol;
y_OL{k_N,k_e} = y_ol;
x_OL{k_N,k_e} = x_ol;

% normalization factors
u_fac = max(abs(u_ol),[],2);
y_fac = max(abs(y_ol),[],2);

% normalize data
u_ol = u_ol./u_fac;
y_ol = y_ol./y_fac;
r    = ref./y_fac;

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
% initialize data structures for controllers
% for k_c = 1:num_c
%     data(k_c).uf_k    = cell(CL_sim_steps,1);
%     data(k_c).yfhat_k = data(k_c).uf_k;
% end

% scale weights (since IO data is also scaled)
Qk_n = y_fac.'*Qk*y_fac;
Rk_n = u_fac.'*Rk*u_fac;
dRk_n= u_fac.'*dRk*u_fac.';

% initialize cell arrays
Cz_run = cell(num_c,1); % controllers
u_run  = cell(num_c,1); % -> u_CL
y_run  = cell(num_c,1); % -> y_CL
x_run  = cell(num_c,1); % -> y_CL

du = mvnrnd(zeros(nu,1),Rdu,CL_sim_steps).';
du_CL{k_N,k_e} = du;

% 1) DeePC with IV
Cz_run{1} =    DeePC(u_ol,y_ol,p,f,N_OL,Qk_n,Rk_n,dRk_n,useAnalytic=true);

% 2) CL-DeePC with IV
Cz_run{2} = CL_DeePC(u_ol,y_ol,p,f,N_CL,Qk_n,Rk_n,dRk_n,useAnalytic=true);

for k_c = 1:num_c
    
    Cz_kc = Cz_run{k_c};
    cost = 0;
    x_k = x0_CL;
    u_last = u_ol(:,end);
    r_all = r;
    e_all = e;
    plant2= plant;

    % initialize u, y, x
    u_run = nan(nu,num_steps); u_run(:,1:Nbar)   =  u_ol;
    y_run = nan(nu,num_steps); y_run(:,1:Nbar)   =  y_ol;
    x_run = nan(nx,num_steps); x_run(:,1:Nbar+1) = [x_ol x_k];

    % set counters
    k1 = Nbar + 1;
    k3 = 1;
    k4 = k3+f-1;
    
    % first computed input
    [uf_k,~] = Cz_kc.solve(rf=r_all(:,k3:k4));
    u_k = uf_k(:,1)+du(:,1);
    u_run(:,k1) = u_k;
    
    % step plant
    e_k = e(:,k1);
    x_k = x_run(:,k1);
    y_k    = plant2.C*x_k + plant2.D*[u_k; e_k];
    x_next = plant2.A*x_k + plant2.B*[u_k; e_k];
    y_run(:,k1) = y_k;

    % update cost
    er_k = y_k-r(:,k3);
    du_k = u_k-u_last;
    cost = cost + er_k.'*Qk_n*er_k + du_k.'*dRk_n*du_k + u_k.'*Rk_n*u_k;
    
    for k2 = Nbar+2:num_steps
        k3 = k3 + 1;
        k4 = k4 + 1;
        
        % previous input
        u_last = u_k;

        % current state
        x_k = x_next;
        x_run(:,k2) = x_k;

        % compute input
        [uf_k,~] = Cz_kc.step(u_k,y_k, rf=r_all(:,k3:k4)); 
        u_k = uf_k(:,1)+du(:,k3);
        u_run(:,k2) = u_k;
        
        % step plant
        e_k    = e(:,k2);
        y_k    = plant2.C*x_k + plant2.D*[u_k; e_k];
        x_next = plant2.A*x_k + plant2.B*[u_k; e_k];
        y_run(:,k2) = y_k;
    
        % update cost
        er_k = y_k-r(:,k3);
        du_k = u_k-u_last;
        cost = cost + er_k.'*Qk_n*er_k + du_k.'*dRk_n*du_k + u_k.'*Rk_n*u_k;
        
        % show progression
        if rem(k3,100)==0
            clc;
            round(k3/CL_sim_steps*100,2)
        end
    end

    % saving data
    u_CL{k_N,k_e,k_c} = u_fac.*u_run(:,end-CL_sim_steps+1:end);
    y_CL{k_N,k_e,k_c} = y_fac.*y_run(:,end-CL_sim_steps+1:end);
    x_CL{k_N,k_e,k_c} =        x_run(:,end-CL_sim_steps+1:end);
    Cz{k_N,k_e,k_c}   = Cz_kc;
    Cost{k_N,k_e,k_c} = cost;
    
end

end
toc
end
%%
plot_all(u_OL,y_OL,u_CL,y_CL,ref,Cz,1,k_N,1,Label,Color)

%%
function plot_all(u_OL,y_OL,u_CL,y_CL,ref,Cz,fignum,k_N,k_e,Label,Color)
    
    figure(fignum);
    clf;
    
    num_c = size(Cz,3);
    Nbar = Cz{k_N,1,1}.Nbar;

    u_ol = u_OL{k_N,k_e};
    y_ol = y_OL{k_N,k_e};

    % subplot with outputs & reference
    ax1 = subplot(2,1,1);
    xline(ax1,2*Nbar+0.5,'k--','HandleVisibility','off'); hold on;
    plot(ax1,Nbar+1:Nbar+length(ref),ref,'--',...
        'DisplayName','reference',...
        'color', [.5 .5 .5], ... grey
        'linewidth', 1.5);%      thicker line
    ylabel(ax1,'$y_k$','interpreter','latex');
    grid on
    
    % subplot with inputs
    ax2 = subplot(2,1,2);
    xline(ax2,2*Nbar+0.5,'k--'); hold on;
    yline(ax2,-15,'r--');
    yline(ax2, 15,'r--');
%     ylim(ax2,[-16,16]);
    ylabel(ax2,'$u_k$','interpreter','latex')
    grid on
    
    lineref1 = cell(num_c,1);
    lineref2 = cell(num_c,1);
    Styles = struct('LineStyle',cell(num_c,1),'Color',cell(num_c,1),'LineWidth',cell(num_c,1));
    for k_c = num_c:-1:1
        
        u_cl = u_CL{k_N,k_e,k_c};
        y_cl = y_CL{k_N,k_e,k_c};
    
        u = [u_ol,u_cl];
        y = [y_ol,y_cl];
                
        lineref1{k_c} = plot(ax1,1:length(y),y,'DisplayName', Label{k_c});
        lineref1{k_c}.LineWidth = 1;
        lineref1{k_c}.Color = Color{k_c};
        Styles(k_c).LineStyle = lineref1{k_c}.LineStyle;
        Styles(k_c).Color     = lineref1{k_c}.Color;
        Styles(k_c).LineWidth = lineref1{k_c}.LineWidth;
        % for k3 = 1:k2
        %     plot(ax1,OL_steps+k3:OL_steps+f+k3-1,cont_struct.yfhat_k{k3})
        % end

        styles2use = namedargs2cell(Styles(k_c));
        lineref2{k_c} = plot(ax2,1:length(u),u,styles2use{:});
        % for k3 = 1:k2
        %     plot(ax2,OL_steps+k3:OL_steps+f+k3-1,cont_struct.uf_k{k3})
        % end
    end
    legend(ax1);

%     ax3 = subplot(3,1,3);
%     plot(ax3,1:length(e),e)
%     xline(ax3,OL_steps+0.5,'k--');
%     ylabel('$e_k$','interpreter','latex')
%     xlabel('samples','interpreter','latex')
%     grid on
    
    linkaxes([ax1 ax2],'x')
    xlim([Nbar,length(u)])
end