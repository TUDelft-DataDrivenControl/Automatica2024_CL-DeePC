function loop_var(x0,N_OL,N_CL,p,f,k_var,k_e,plant,Ru,Re,ny,nu,nx,num_steps,Nbar,ref,Qk,Rk,dRk,num_c,Rdu,CL_sim_steps,dir_name,seed_num,Obsv_f,Lu_act,Ly_act,Gu_act)

%% user defined constraints
y_max = 1000;
u_max = 15;
du_max = 3.75;

%% make used stochastic data
rng(seed_num);
% -> save now s.t. if there is an error this data is available for debugging
OL_sim_steps = num_steps-CL_sim_steps;
e     = mvnrnd(zeros(ny,1),Re, num_steps).';    % noise realization
u_OL  = mvnrnd(zeros(nu,1),Ru, OL_sim_steps).'; % OL input
du_CL = mvnrnd(zeros(nu,1),Rdu,CL_sim_steps).'; % CL input

% maximize input disturbance so as not to cause infeasibility
du_CL(du_CL<-du_max) = -du_max; % shouldn't affect many instances -> otherwise changes effective Rdu
du_CL(du_CL>du_max)  =  du_max;

% saving data
name_kvar = inputname(6);
save_str = fullfile(dir_name, sprintf('%s_%d_ke_%d_ks_%d.mat',name_kvar([1 3]),k_var,k_e,seed_num));
save(save_str,'e','u_OL','du_CL');

% splitting e into OL & CL parts
e_OL = e(:,1:OL_sim_steps);
e_CL = e(:,OL_sim_steps:end);

%% initial open loop simulation

% simulate system
[y_OL,~,x_OL] = lsim(plant,[u_OL;e_OL],[],x0);
y_OL = y_OL.'; x_OL = x_OL.';
x0_CL = plant.A*x_OL(:,end) + plant.B*[u_OL(:,end); e_OL(:,end)];

% select used OL data
u_ol = u_OL(:,end-Nbar+1:end);
y_ol = y_OL(:,end-Nbar+1:end);

%% closed-loop operation

% step-simulate system
simulate_step = @(u_k,e_k,x_k) deal(plant.C*x_k + plant.D*[u_k;e_k],... [y_k,
                                    plant.A*x_k + plant.B*[u_k;e_k]);  % x_next]
stage_cost    = @(u,y_k) y_k.'*Qk*y_k + u(:,2).'*Rk*u(:,2) + (u(:,2)-u(:,1)).'*dRk*(u(:,2)-u(:,1));

% user defined constraints
con = struct();
con.y_max  =  y_max;
con.y_min  = -y_max;
con.u_max  =  u_max;
con.u_min  = -u_max;
con.du_max = du_max;

% initialize data arrays
Cz = cell(num_c,1);
[u_CL,y_CL,x_CL,Cost,eLu,eLy,eGu,eObX,stat] = deal(Cz);

% 1) DeePC with IV
Cz{1} =    DeePC(u_ol,y_ol,p,f,N_OL,Qk,Rk,dRk,constr=con);

% 2) CL-DeePC with IV
Cz{2} = CL_DeePC(u_ol,y_ol,p,f,N_CL,Qk,Rk,dRk,constr=con);

for k_c = 1:num_c

    % initialize data for run with controller
    x_CLr = nan(nx,CL_sim_steps+1); x_CLr(:,1) = x0_CL;
    [u_CLr,y_CLr,cost,eLu_r,eLy_r,eGu_r,eObX_r,stat_r] = deal( nan(nu,CL_sim_steps) );
    
    % get first CL input
    [uf_k,~,stat_r(1)] = Cz{k_c}.solve(rf=ref(:,1:f));
    % add disturbance input
    u_CLr(:,1) = uf_k(:,1)+du_CL(:,1);
    % simulate step
    [y_CLr(:,1),x_CLr(:,2)] = simulate_step(u_CLr(:,1),e_CL(:,1),x_CLr(:,1));
    
    % analysis
    cost(1) = stage_cost([u_OL(:,end) u_CLr(:,1)],y_CLr(:,1)); % stage cost
    [eObX_r(1),eLu_r(1), eLy_r(1), eGu_r(1)] = error_ID(Cz{k_c},x_CLr(:,1),Obsv_f,Lu_act,Ly_act,Gu_act);      % ID errors


    for k = 2:CL_sim_steps
        % get input
        try
            [uf_k,~,stat_r(k)] = Cz{k_c}.step(u_CLr(:,k-1), y_CLr(:,k-1), rf=ref(:,k:k+f-1));
        catch Error
            disp(['k_var =',num2str(k_var),'; k_e = ',num2str(k_e),' k1 = ',num2str(k)]);
            error(Error.message)
        end
        % add disturbance input
        u_CLr(:,k) = uf_k(:,1)+du_CL(:,k);
        % simulate step
        [y_CLr(:,k),x_CLr(:,k+1)] = simulate_step(u_CLr(:,k),e_CL(:,k),x_CLr(:,k));
        % calculate cost
        cost(k) = stage_cost(u_CLr(:,k-1:k),y_CLr(:,k));      % stage cost
        [eObX_r(k),eLu_r(k), eLy_r(k), eGu_r(k)] = error_ID(Cz{k_c},x_CLr(:,k),Obsv_f,Lu_act,Ly_act,Gu_act); % ID errors
%         if k==2 || k/CL_sim_steps*100>=k2
%             if k==2; k2 = 10; else k2 = k2 + 5; end
%             clc
%             disp(strcat('Progress: k_e = ',sprintf('%d',k_e),', Controller #',num2str(k_c),' , ',num2str(k/CL_sim_steps*100),'%'));
%         end
    end

    % saving CL data from run with controller
    u_CL{k_c} = u_CLr;
    y_CL{k_c} = y_CLr;
    x_CL{k_c} = x_CLr;
    Cost{k_c} = cost;
    eLu{k_c}  = eLu_r;
    eLy{k_c}  = eLy_r;
    eGu{k_c}  = eGu_r;
    eObX{k_c} = eObX_r;
    stat{k_c} = stat_r;

end % end for k_c

save(save_str,'y_OL','x_OL','u_CL','y_CL','x_CL','Cost','eLu','eLy','eGu','eObX','stat','-append')
end

%% Helper functions

% limiter of u: rate of change & min/max
function u_k = limiter(u_k,u_last,u_max,du_max)
    if any(abs(u_k) > u_max) || any(abs(u_k-u_last) > du_max)
        u_max = min([u_max,u_last+du_max],[],2);
        u_min = max([-u_max,u_last-du_max],[],2);
        u_k = min([max([u_k,u_min],[],2),u_max],[],2);
    end
end

% ID error - from Lu, Ly, Gu
function [eObX,eLu,eLy,eGu] = error_ID(Cz,x_k,Obsv_f,Lu_act,Ly_act,Gu_act)
    eObX = norm(Cz.Prob.Lu*Cz.up(:)+Cz.Prob.Ly*Cz.yp(:)-Obsv_f*x_k,'Inf');
    eLu = norm(Cz.Prob.Lu-Lu_act,'fro')^2;
    eLy = norm(Cz.Prob.Ly-Ly_act,'fro')^2;
    eGu = norm(Cz.Prob.Gu-Gu_act,'fro')^2;
end