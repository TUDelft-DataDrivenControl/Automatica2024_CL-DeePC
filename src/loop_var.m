function loop_var(x0,N_OL,N_CL,p,f,k_var,k_e,plant,Ru,Re,ny,nu,nx,num_steps,Nbar,ref,Qk,Rk,dRk,num_c,Rdu,CL_sim_steps,temp_str)

%% initial open loop simulation

% noise realization
e  = mvnrnd(zeros(ny,1),Re,num_steps).';
u_ol = mvnrnd(zeros(nu,1),Ru,Nbar).';

% simulate system
[y_ol,~,x_ol] = lsim(plant,[u_ol;e(:,1:Nbar)],[],x0);
y_ol = y_ol.'; x_ol = x_ol.';
x0_CL = plant.A*x_ol(:,end) + plant.B*[u_ol(:,end); e(:,Nbar)];

% saving data - OL
u_OL = u_ol;
y_OL = y_ol;
x_OL = x_ol;

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

% input disturbance
du = mvnrnd(zeros(nu,1),Rdu,CL_sim_steps).';
du_CL = du;
du = du./u_fac;

du_max = 3.75;
con = struct('Opti',cell(1,2));
for k_con = 1:num_c
con(k_con).Opti = casadi.Opti();
con(k_con).uf   = con(k_con).Opti.variable(nu,f);
con(k_con).u0   = con(k_con).Opti.parameter(nu,1);
con(k_con).expr = {con(k_con).uf <=  15./u_fac;
            con(k_con).uf >= -15./u_fac;
            con(k_con).u0 - con(k_con).uf(:,1) <=  du_max./u_fac;
            con(k_con).u0 - con(k_con).uf(:,1) >= -du_max./u_fac;
            con(k_con).uf(:,1:end-1)-con(k_con).uf(:,2:end) <=  du_max./u_fac;
            con(k_con).uf(:,1:end-1)-con(k_con).uf(:,2:end) >= -du_max./u_fac};
end

% 1) DeePC with IV
Cz{1} =    DeePC(u_ol,y_ol,p,f,N_OL,Qk_n,Rk_n,dRk_n,constr=con(1),UseOptimizer=true);

% 2) CL-DeePC with IV
Cz{2} = CL_DeePC(u_ol,y_ol,p,f,N_CL,Qk_n,Rk_n,dRk_n,constr=con(2),UseOptimizer=true);

for k_c = 1:num_c
    cost = 0;
    x_k = x0_CL;
    u_last = u_ol(:,end);

    % initialize u, y, x
    u_run = nan(nu,num_steps); u_run(:,1:Nbar)   =  u_ol;
    y_run = nan(nu,num_steps); y_run(:,1:Nbar)   =  y_ol;
    x_run = nan(nx,num_steps); x_run(:,1:Nbar+1) = [x_ol x_k];

    % set counters
    k1 = Nbar + 1;
    k3 = 1;
    k4 = k3+f-1;
    
    % first computed input
    [uf_k,~] = Cz{k_c}.solve(rf=r(:,k3:k4));
    u_k = uf_k(:,1)+du(:,1);
    u_run(:,k1) = u_k;
    
    % step plant
    e_k = e(:,k1);
    x_k = x_run(:,k1);
    y_k    = plant.C*x_k + plant.D*[u_k; e_k];
    x_next = plant.A*x_k + plant.B*[u_k; e_k];
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
        [uf_k,~] = Cz{k_c}.step(u_k, y_k, rf=r(:,k3:k4)); 
        u_k = uf_k(:,1)+du(:,k3);
        u_run(:,k2) = u_k;
        
        % step plant
        e_k    = e(:,k2);
        y_k    = plant.C*x_k + plant.D*[u_k; e_k];
        x_next = plant.A*x_k + plant.B*[u_k; e_k];
        y_run(:,k2) = y_k;
    
        % update cost
        er_k = y_k-r(:,k3);
        du_k = u_k-u_last;
        cost = cost + er_k.'*Qk_n*er_k + du_k.'*dRk_n*du_k + u_k.'*Rk_n*u_k;
    end
    % saving data - CL
    u_CL{k_c} = u_fac.*u_run(:,end-CL_sim_steps+1:end); %{k_var,k_e,k_c}
    y_CL{k_c} = y_fac.*y_run(:,end-CL_sim_steps+1:end);
    x_CL{k_c} =        x_run(:,end-CL_sim_steps+1:end);
    Cost{k_c} = cost;

end % end for k_c

save(strcat('../data/temp/',temp_str,num2str(k_var),'_ke_',num2str(k_e),'.mat'),'u_OL','y_OL','x_OL','e','du_CL','u_CL','y_CL','x_CL','Cost')