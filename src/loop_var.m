function loop_var(x0,N_OL,N_CL,p,f,k_var,k_e,plant,Ru,Re,ny,nu,nx,num_steps,Nbar,ref,Qk,Rk,dRk,num_c,Rdu,CL_sim_steps,temp_str,seed_num)

%% make used stochastic data
rng(seed_num);
% -> save now s.t. if there is an error this data is available for debugging
OL_sim_steps = num_steps-CL_sim_steps;
e  = mvnrnd(zeros(ny,1),Re,num_steps).';        % noise realization
u_OL = mvnrnd(zeros(nu,1),Ru,OL_sim_steps).';   % OL input
du_CL = mvnrnd(zeros(nu,1),Rdu,CL_sim_steps).'; % CL input

% saving data
save_str = strcat('../data/temp/',temp_str,num2str(k_var),'_ke_',num2str(k_e),'.mat');
save(save_str,'e','u_OL','du_CL');

%% initial open loop simulation

% simulate system
[y_OL,~,x_OL] = lsim(plant,[u_OL;e(:,1:OL_sim_steps)],[],x0);
y_OL = y_OL.'; x_OL = x_OL.';
x0_CL = plant.A*x_OL(:,end) + plant.B*[u_OL(:,end); e(:,OL_sim_steps)];

% select used OL data
u_ol = u_OL(:,end-Nbar+1:end);
y_ol = y_OL(:,end-Nbar+1:end);
x_ol = x_OL(:,end-Nbar+1:end);

%% closed-loop operation

% user defined constraints
y_max = 1000;
u_max = 15;
du_max = 3.75;
con = struct();
con.y_max  =  y_max;
con.y_min  = -y_max;
con.u_max  =  u_max;
con.u_min  = -u_max;
con.du_max = du_max;

% initialize arrays
Cz = cell(num_c,1);
u_CL = Cz;
y_CL = Cz;
x_CL = Cz;
Cost = Cz;

% 1) DeePC with IV
Cz{1} =    DeePC(u_ol,y_ol,p,f,N_OL,Qk,Rk,dRk,constr=con);

% 2) CL-DeePC with IV
Cz{2} = CL_DeePC(u_ol,y_ol,p,f,N_CL,Qk,Rk,dRk,constr=con);

for k_c = 1:num_c
    cost = 0;
    x_k = x0_CL;
    u_last = u_ol(:,end);

    % initialize u, y, x
    u_run = nan(nu,Nbar+CL_sim_steps); u_run(:,1:Nbar)   =  u_ol;
    y_run = nan(nu,Nbar+CL_sim_steps); y_run(:,1:Nbar)   =  y_ol;
    x_run = nan(nx,Nbar+CL_sim_steps); x_run(:,1:Nbar+1) = [x_ol x_k];

    % set counters
    k1 = Nbar + 1;
    k2 = k1+f-1;
    k3 = 1;
    
    % first computed input
    [uf_k,~] = Cz{k_c}.solve(rf=ref(:,k1:k2));
    u_k = limiter(uf_k(:,1),u_last,u_max,du_max);
    u_k = u_k+du_CL(:,1);
    u_run(:,k1) = u_k;
    
    % step plant
    e_k = e(:,k1);
    x_k = x_run(:,k1);
    y_k    = plant.C*x_k + plant.D*[u_k; e_k];
    x_next = plant.A*x_k + plant.B*[u_k; e_k];
    y_run(:,k1) = y_k;

    % update cost
    er_k = y_k-ref(:,k3);
    du_k = u_k-u_last;
    cost = cost + er_k.'*Qk*er_k + du_k.'*dRk*du_k + u_k.'*Rk*u_k;
    
    for k1 = Nbar+2:Nbar+CL_sim_steps
        k2 = k2 + 1;
        k3 = k3 + 1;
        
        % previous input
        u_last = u_k;

        % current state
        x_k = x_next;
        x_run(:,k1) = x_k;

        % compute input
        try
            [uf_k,~] = Cz{k_c}.step(u_k, y_k, rf=ref(:,k1:k2));
        catch Error
            disp(['k_var =',num2str(k_var),'; k_e = ',num2str(k_e),' k1 = ',num2str(k3)]);
            error(Error.message)
        end
        u_k = limiter(uf_k(:,1),u_last,u_max,du_max);
        u_k = u_k+du_CL(:,k3);
        u_run(:,k1) = u_k;
        
        % step plant
        e_k    = e(:,k1);
        y_k    = plant.C*x_k + plant.D*[u_k; e_k];
        x_next = plant.A*x_k + plant.B*[u_k; e_k];
        y_run(:,k1) = y_k;
    
        % update cost
        er_k = y_k-ref(:,k3);
        du_k = u_k-u_last;
        cost = cost + er_k.'*Qk*er_k + du_k.'*dRk*du_k + u_k.'*Rk*u_k;

%         if rem(k3,100) == 0
%             disp(strcat('ke',num2str(k_e),' progress: ',num2str(k3/CL_sim_steps*100),'%'))
%         end
    end
    % saving data - CL
    u_CL{k_c} = u_run(:,end-CL_sim_steps+1:end);
    y_CL{k_c} = y_run(:,end-CL_sim_steps+1:end);
    x_CL{k_c} = x_run(:,end-CL_sim_steps+1:end);
    Cost{k_c} = cost;

end % end for k_c

save(save_str,'y_OL','x_OL','u_CL','y_CL','x_CL','Cost','-append')
end

function u_k = limiter(u_k,u_last,u_max,du_max)
    if any(abs(u_k) > u_max) || any(abs(u_k-u_last) > du_max)
        u_max = min([u_max,u_last+du_max],[],2);
        u_min = max([-u_max,u_last-du_max],[],2);
        u_k = min([max([u_k,u_min],[],2),u_max],[],2);
    end
end