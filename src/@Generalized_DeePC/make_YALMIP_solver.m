function make_YALMIP_solver(obj,usr_con,expr_flag)
% construct solve method using YALMIP

str2 = {'u_min','u_max','y_min','y_max','du_max','dy_max'};

%% process usr_con.expr if it is a cell array
if contains(expr_flag,'cell')
    % obtain con_LHS & con_gleq from usr_con.expr
    if ~isvector(usr_con.expr) && ismatrix(usr_con.expr)
        con_LHS  = usr_con.expr(:,1);
        con_gleq = usr_con.expr(:,2);
    elseif isvector(usr_con.expr)
        con_LHS  = usr_con.expr(1:2:end); con_LHS = con_LHS(:);
        con_gleq = usr_con.expr(2:2:end); con_gleq= con_gleq(:);
    end

    usr_con.expr = [];
    for k = 1:numel(con_LHS)
        if strcmp(con_gleq{k},'==')
            usr_con.expr = [usr_con.expr;con_LHS{k}==0];
        elseif strcmp(con_gleq{k},'<=') || strcmp(con_gleq{k},'=<')
            usr_con.expr = [usr_con.expr;con_LHS{k}<=0];
        elseif strcmp(con_gleq{k},'>=') || strcmp(con_gleq{k},'=>')
            usr_con.expr = [usr_con.expr;con_LHS{k}>=0];
        else
            error('Constraint specified incorrectly.')
        end
    end
end

obj.Prob.con_usr = usr_con.expr;

%% process usr_con str2 fields
for k=1:length(str2)
    fn = str2{k};
    if isfield(usr_con,fn)
        switch k
            case 1 % u_min
                usr_con.expr = [usr_con.expr; obj.Prob.uf_(:) >= repmat(usr_con.u_min,obj.f,1)];
            case 2 % u_max
                usr_con.expr = [usr_con.expr; obj.Prob.uf_(:) <= repmat(usr_con.u_max,obj.f,1)];
            case 3 % y_min
                usr_con.expr = [usr_con.expr; obj.Prob.yf_(:) >= repmat(usr_con.y_min,obj.f,1)];
            case 4 % y_max
                usr_con.expr = [usr_con.expr; obj.Prob.yf_(:) <= repmat(usr_con.y_max,obj.f,1)];
            case 5 % du_max
                du_ = [obj.Prob.uf_(:,1)-obj.Prob.up_(:,end) obj.Prob.uf_(:,2:end)-obj.Prob.uf_(:,1:end-1)];
                du_max = repmat(usr_con.du_max,obj.f,1);
                usr_con.expr = [usr_con.expr; du_(:) >= -du_max; du_(:) <= du_max];
            case 6 % dy_max
                dy_ = [obj.Prob.yf_(:,1)-obj.Prob.yp_(:,end) obj.Prob.yf_(:,2:end)-obj.Prob.yf_(:,1:end-1)];
                dy_max = repmat(usr_con.dy_max,obj.f,1);
                usr_con.expr = [usr_con.expr; dy_(:) >= -dy_max; dy_(:) <= dy_max];
        end
    end
end

%% make constraints for dynamics
if obj.options.ExplicitPredictor
    obj.make_con_dyn_ExplicitPredictor();
else
    obj.make_con_dyn_4Optimization();
end

%% Make solver method

% solver options
opts = obj.options.opts;
if isempty(opts)
    obj.Prob.sdp_opts = namedargs2cell(obj.Prob.sdp_opts);
    obj.Prob.sdp_opts = sdpsettings(obj.Prob.sdp_opts{:});
else
    if ~isfield(opts,'solver')
        opts.solver = obj.Prob.sdp_opts.solver;
    end
    opts = namedargs2cell(opts);
    obj.Prob.sdp_opts = sdpsettings(opts);
end

% make solver
if obj.options.UseOptimizer
    constraints = [obj.Prob.con_dyn;...
                   obj.Prob.con_usr];
    if obj.options.ExplicitPredictor
        obj.Prob.Optimizer = ...
        optimizer(constraints,obj.Prob.cost,obj.Prob.sdp_opts,...
                 {obj.Prob.Lu_,obj.Prob.Ly_,obj.Prob.Gu_,obj.Prob.up_,obj.Prob.yp_,obj.Prob.rf_},... Parameters
                 {obj.Prob.uf_,obj.Prob.yf_});  % Outputs
    else
        obj.Prob.Optimizer = ...
        optimizer(constraints,obj.Prob.cost,obj.Prob.sdp_opts,...
                {obj.Prob.LHS_,obj.Prob.up_,obj.Prob.yp_,obj.Prob.rf_},... Parameters
                {obj.Prob.uf_,obj.Prob.yf_});  % Outputs
    end
    obj.solve = @obj.optimizer_solve;
else
    obj.solve = @obj.optimize_solve;
end

