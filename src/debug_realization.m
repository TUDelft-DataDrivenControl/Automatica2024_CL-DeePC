close all;
clc;

%% Simulation settings
model_Favoreel1999 % loads model from Favoreel 1999 - original SPC paper

rho_max = max(abs(eig(A-K*C)),[],'all');
cond1_fac = -1/(2*log(rho_max)); % cond1_fac*log(N) < p

Nmax = 10^3;                      % maximum number of columns
pmin = ceil(cond1_fac*log(Nmax)); % such that above always true

% controller settings
p = 20; % > pmin with Nmax
f = 20;
Nmin = p*(nu+ny)+f*nu; % minimum N is determined by regular DeePC
Qk = 10;
Rk = 0.1;
dRk= 0;

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
Ru = 1*eye(nu);    % OL input
Rdu= Ru/4;         % CL input disturbance

% OL-sim initial state
x0 = zeros(nx,1);

% number of CL simulation steps
OL_sim_steps = 1200; %>=Nbar_max
CL_sim_steps = 1800;
num_steps = OL_sim_steps + CL_sim_steps;  % total simulation length

% define reference
% ref = nan(ny,CL_sim_steps+f-1); % +f-1 needed for simulation end
% ref = (-square((0:size(ref,2)-1)*2*pi/(200)))*50+1*50;
ref        = 50*[-sign(sin(2*pi/2*0.01*(0:CL_sim_steps-1))) ones(1,f-1)]+500;

%%
k_N = 1;
N_OL = N_all_OL(k_N); N_OL = 500;
N_CL = N_all_CL(k_N); N_CL = N_OL+f-1;
Nbar = Nbar_all(k_N); Nbar = N_OL+p+f-1;
%%
unoise = mvnrnd(zeros(nu,1),Ru,num_steps+1).';
e      = mvnrnd(zeros(ny,1),Re,num_steps+1).';
%% initial open loop simulation

% OL_sim_steps = num_steps-CL_sim_steps;
u_OL = unoise(:,1:OL_sim_steps); du_CL = 0.5*unoise(:,OL_sim_steps+1:end);
% simulate system
[y_ol,~,x_ol] = lsim(plant,[u_OL;e(:,1:OL_sim_steps)],[],x0);
y_ol = y_ol.'; x_ol = x_ol.';
x0_CL = plant.A*x_ol(:,end) + plant.B*[u_OL(:,end); e(:,OL_sim_steps)];

% differentiate between normalized & original OL data
u_ol = u_OL;
y_OL = y_ol;
x_OL = x_ol;

% normalization factors
u_fac = 1;%max([max(abs(u_ol),[],2),15]);
y_fac = 1;%max([max(abs(y_ol),[],2),100]);

% normalize data
u_ol = u_ol./u_fac;
y_ol = y_ol./y_fac;
r    = ref./y_fac;

% truncate data
u2_ol = u_ol(:,end-Nbar+1:end);
y2_ol = y_ol(:,end-Nbar+1:end);
x2_ol = x_ol(:,end-Nbar+1:end);

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
% plant.B(:,1:nu)     =             plant.B(:,1:nu)*diag(u_fac);
% plant.C             = diag(y_fac)\plant.C;
% plant.D(:,1:nu)     = diag(y_fac)\plant.D(:,1:nu)*diag(u_fac);
% plant.D(:,nu+1:end) = diag(y_fac)\plant.D(:,nu+1:end);

%% closed-loop operation
% scale weights (since IO data is also scaled)
Qk_n = y_fac.'*Qk*y_fac;
Rk_n = u_fac.'*Rk*u_fac;
dRk_n= u_fac.'*dRk*u_fac.';

% initialize arrays
Cz = cell(num_c,1);
u_CL = Cz;
y_CL = Cz;
x_CL = Cz;
Cost = Cz;

% user defined constraints
u_max = 15;     u_max =  u_max./u_fac;
du_max = 3.75; du_max = du_max./u_fac;

qpoases_opts = struct('solver','qpoases','options',struct('printLevel','none'));
nlpsol_opts  = struct('solver','nlpsol','options',struct('nlpsol','ipopt','nlpsol_options',struct('ipopt',struct('print_level',0,'warm_start_init_point','yes','nlp_scaling_method','none'),'print_time',0)));
ipopt_opts   = struct('solver','ipopt','options',struct('ipopt',struct('print_level',0,'warm_start_init_point','no','nlp_scaling_method','none','fixed_variable_treatment','make_constraint'),'print_time',0));

con = struct();
con.u_min  = -u_max;
con.u_max  =  u_max;
con.y_max  =  1000./y_fac;
con.y_min  = -1000./y_fac;
con.du_max = du_max;
Cz{1} =    DeePC(u2_ol,y2_ol,p,f,N_OL,Qk_n,Rk_n,dRk_n,constr=con);%opts=nlpsol_opts);
Cz{2} = CL_DeePC(u2_ol,y2_ol,p,f,N_CL,Qk_n,Rk_n,dRk_n,constr=con);%opts=nlpsol_opts);

% input disturbance
du = du_CL./u_fac;

for k_c = 1:num_c
    cost1 = 0;
    x_k = x0_CL;
    u_last = u_ol(:,end);

    % initialize u, y, x
    u_run = nan(nu,OL_sim_steps+CL_sim_steps); u_run(:,1:OL_sim_steps)   =  u_ol;
    y_run = nan(nu,OL_sim_steps+CL_sim_steps); y_run(:,1:OL_sim_steps)   =  y_ol;
    x_run = nan(nx,OL_sim_steps+CL_sim_steps); x_run(:,1:OL_sim_steps+1) = [x_ol x_k];

    % set counters
    k1 = OL_sim_steps + 1;
    k2 = 1;
    k3 = k2+f-1;
    
    % first computed input
    [uf_k,~] = Cz{k_c}.solve(rf=r(:,k1:k1+f-1));
    u_k = limiter(uf_k(:,1),u_last,u_max,du_max);
    u_k = u_k+du(:,1);
    u_run(:,k1) = u_k;
    
    % step plant
    e_k = e(:,k1);
    x_k = x_run(:,k1);
    y_k    = plant.C*x_k + plant.D*[u_k; e_k];
    x_next = plant.A*x_k + plant.B*[u_k; e_k];
    y_run(:,k1) = y_k;

    % update cost
    er_k = y_k-r(:,k2);
    du_k = u_k-u_last;
    cost1 = cost1 + er_k.'*Qk_n*er_k + du_k.'*dRk_n*du_k + u_k.'*Rk_n*u_k;
    
    for k1 = OL_sim_steps+2:OL_sim_steps+CL_sim_steps
        k2 = k2 + 1;
        k3 = k3 + 1;
        
        % previous input
        u_last = u_k;

        % current state
        x_k = x_next;
        x_run(:,k1) = x_k;

        % compute input
        try
            [uf_k,~] = Cz{k_c}.step(u_k, y_k, rf=r(:,k2:k3));
        catch Error
            figure()
            ax1 = subplot(2,1,1);
            plot(y_fac.*y_run);
            hold on; grid on;
            plot(OL_sim_steps+1:OL_sim_steps+length(ref),ref);
            ax2 = subplot(2,1,2);
            plot(u_fac.*u_run);
            grid on
            linkaxes([ax1 ax2],'x')

            disp(['debug k2 = ',num2str(k2)]);
            error(Error.message)
        end
        u_k = limiter(uf_k(:,1),u_last,u_max,du_max);
        u_k = u_k+du(:,k2);
        u_run(:,k1) = u_k;
        
        % step plant
        e_k    = e(:,k1);
        y_k    = plant.C*x_k + plant.D*[u_k; e_k];
        x_next = plant.A*x_k + plant.B*[u_k; e_k];
        y_run(:,k1) = y_k;
    
        % update cost
        er_k = y_k-r(:,k2);
        du_k = u_k-u_last;
        cost1 = cost1 + er_k.'*Qk_n*er_k + du_k.'*dRk_n*du_k + u_k.'*Rk_n*u_k;

        if rem(k2,100) == 0
            disp(strcat('progress: ',num2str(k2/CL_sim_steps*100),'%'))
        end
    end
    % saving data - CL
    u_CL{k_c} = u_fac.*u_run(:,end-CL_sim_steps+1:end);
    y_CL{k_c} = y_fac.*y_run(:,end-CL_sim_steps+1:end);
    x_CL{k_c} =        x_run(:,end-CL_sim_steps+1:end);
    Cost{k_c} = cost1;
      
    figure()
    ax1 = subplot(2,1,1);
    plot(y_fac.*y_run);
    hold on; grid on;
    plot(OL_sim_steps+1:OL_sim_steps+length(ref),ref);
    ax2 = subplot(2,1,2);
    plot(u_fac.*u_run);
    grid on
    linkaxes([ax1 ax2],'x')

end % end for k_c

function u_k = limiter(u_k,u_last,u_max,du_max)
    if any(abs(u_k) > u_max) || any(abs(u_k-u_last) > du_max)
        u_max = min([u_max,u_last+du_max],[],2);
        u_min = max([-u_max,u_last-du_max],[],2);
        u_k = min([max([u_k,u_min],[],2),u_max],[],2);
    end
end
