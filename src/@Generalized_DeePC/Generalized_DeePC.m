classdef Generalized_DeePC < handle
    %CL_DEEPC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        p
        f
        N
        Nbar
        nu
        ny

        options = struct('use_IV',[],'adaptive',[],'ExplicitPredictor',false)
        solve            % method has yalmip optimizer/optimize variant
        step_data_update % method has adaptive/non-adaptive variant

        % data
        Upast
        Ypast
        rf

        % optimization variables & functions
        Prob = struct('cost', [],'yf_',[],'uf_',[],'rf_',[],...
                      'yp_',  [],'up_',[],... % 
                      'G_',   [],'Optimizer',[],...
                      'Hf_',[],'LHS_',[],...
                      'con_usr',[],'con_dyn',[],...
                      'sdp_opts',sdpsettings('solver','mosek','verbose',0),...
                      'cas_opts',struct('solver', 'osqp',...
                                        'options', struct('print_time',0)));
%                                   struct('solver', 'ipopt',...
%                                         'options', struct('print_time',1, ...
%                                                           'ipopt',struct('print_level',0,'nlp_scaling_method','none','warm_start_init_point','no'))));
    end
    properties (Dependent)
        LHS
        up
        yp
    end
    properties(Hidden)
        fid
        pfid
        nGcols
        make_var
        make_par
    end
    
    methods
        function obj = Generalized_DeePC(u,y,p,f,fid,N,Q,R,dR,options,con_user,solve_type)
            arguments
                u   (:,:) double
                y   (:,:) double
                p   (1,1) double {mustBeInteger,mustBePositive}
                f   (1,1) double {mustBeInteger,mustBePositive}
                fid (1,1) double {mustBeInteger,mustBePositive,mustBeLessThanOrEqual(fid,f)}
                N   (1,1) double {mustBeInteger,mustBePositive}
                Q   (:,:) double
                R   (:,:) double
                dR  (:,:) double
                options.use_IV logical   = true
                options.adaptive logical = true
                options.useAnalytic logical = true  % use analytic solution if there are no constraints
                options.ExplicitPredictor logical = true
                con_user.constr struct = struct('expr',[],'u0',[],'uf',[],'y0',[],'yf',[]);
                solve_type.UseOptimizer logical = true
                solve_type.opts = []
            end

            obj.nu = min(size(u)); if size(u,2) < size(u,1); u=u.'; end
            obj.ny = min(size(y)); if size(y,2) < size(y,1); y=y.'; end
            
            obj.p = p;
            obj.f = f;
            
            obj.fid    = fid;
            obj.pfid   = obj.p + obj.fid;
            obj.nGcols = obj.f + 1 - obj.fid;

            % N >= (p+fid)*nu + p*ny <== Guarantee # cols >= # rows
            Nmin = (obj.pfid)*obj.nu + obj.p*obj.ny;
            if N < Nmin
                obj.N = Nmin;
                error('Data matrix must have at least as many columns as rows.')
            else
                obj.N = N;
            end
            obj.Nbar = obj.N+obj.pfid-1;
            
            % parse options
            options_fn = fieldnames(options);
            for k = 1:length(options_fn)
                fn = options_fn{k};
                if ~isempty(options.(fn))
                    obj.options.(fn) = options.(fn);
                end
            end
            
            % initialize data
            validateattributes(u, 'double',{'nrows',obj.nu,'ncols',obj.Nbar})
            validateattributes(y, 'double',{'nrows',obj.ny,'ncols',obj.Nbar})
            if options.adaptive
                obj.Upast = u;
                obj.Ypast = y;
            else
                % 1st part is for constant past data, latter for changing up & yp
                obj.Upast = [u u(:,end-obj.p+1:end)];
                obj.Ypast = [y y(:,end-obj.p+1:end)];
            end
            
            % size cost function weighting matrices
            if size(Q,1)==obj.ny && size(Q,2)==size(Q,1)
                Q = kron(eye(obj.f),Q);
            elseif ~(size(Q,1)==obj.ny*obj.f && size(Q,2)==size(Q,1))
                error('Output error weighting matrix Q of incompatible dimensions.')
            end
            if size(R,1)==obj.nu && size(R,2)==size(R,1)
                R = kron(eye(obj.f),R);
            elseif ~(size(R,1)==obj.nu*obj.f && size(R,2)==size(R,1))
                error('Input weighting matrix R of incompatible dimensions.')
            end
            if size(dR,1)==obj.nu && size(dR,2)==size(dR,1)
                dR = kron(eye(obj.f),dR);
            elseif ~(size(dR,1)==obj.nu*obj.f && size(dR,2)==size(dR,1))
                error('Input weighting matrix dR of incompatible dimensions.')
            end
                % make all weighting matrices symmetric
            [Q,R,dR] = deal((Q+Q.')/2,(R+R.')/2,(dR+dR.')/2);
            [obj.Prob.Q,obj.Prob.R,obj.Prob.dR] = deal(Q,R,dR);

            %==============================================================
            %-------------- make optimization problem ---------------------
            %==============================================================
            if ~isempty(con_user.constr.expr) || ~obj.options.useAnalytic
                % user-defined constraints involved -> use optimization
                obj.Prob.con_usr = con_user.constr.expr;
                usr_con = rmfield(con_user.constr,'expr');
                names = fieldnames(usr_con);
                if isa(usr_con.(names{1}),'sdpvar')
                    obj.options.SolverFramework = 1; % 1 <- use YALMIP
                    obj.make_var = @(dim1,dim2) sdpvar(dim1,dim2,'full');
                    obj.make_par = @(dim1,dim2) sdpvar(dim1,dim2,'full');
                elseif contains(class(usr_con.(names{1})),'casadi')
                    obj.options.SolverFramework = 2; % 2 <- using CasADi
                    if isfield(usr_con,'Opti')
                        obj.Prob.Opti = usr_con.Opti;%.copy();
                        obj.make_var = @(dim1,dim2) obj.Prob.Opti.variable(dim1,dim2);
                        obj.make_par = @(dim1,dim2) obj.Prob.Opti.parameter(dim1,dim2);
                    else
                        error(['Specify the employed instance of the casadi.Opti...' ...
                            ' class as field ".Opti" of the structure containing...' ...
                            ' the user-defined constraints'])
                    end
                else
                    error('Optimization framework not recognized');
                end

                % parse solver options
                if obj.options.SolverFramework == 1 % <- using YALMIP
                    if ~isempty(solve_type.opts) % otherwise use default solver options
                        if ~isfield(solve_type.opts,'solver')
                            solve_type.solver = obj.Prob.sdp_opts.solver;
                        end
                        obj.Prob.sdp_opts = solve_type.opts;
                    end
                else % using CasADi
                    if isempty(solve_type.opts) % use default solver options
                        obj.Prob.Opti.solver(...
                            obj.Prob.cas_opts.solver,...
                            obj.Prob.cas_opts.options);
                    else
                        if ~isfield(solve_type.opts,'solver')
                            solve_type.opts.solver = obj.Prob.cas_opts.solver;
                        end
                        obj.Prob.Opti.solver(solve_type.opts.solver,solve_type.opts.options);
                    end
                end
                
                % parse user-defined constraints.
                names = {'uf','u0','yf','y0'};
                for kk = 1:length(names)
                    name = names{kk};
                    if strcmp(name(1),'u')
                        dim1 = obj.nu;
                    else
                        dim1 = obj.ny;
                    end
                    if isfield(usr_con,name) % -> variable is specified (at least in part) so use
                        user_var = usr_con.(name);
                        if strcmp(name(2),'f') % uf or yf is specified
                            obj.Prob.([name,'_'])=user_var;
                        else % specified u0 or y0 -> last sample of up or yp
                            obj.Prob.([name(1),'p_'])=[obj.make_par(dim1,obj.p-1),user_var];
                        end
                    elseif strcmp(name(2),'f') % define uf or yf parameters
                        obj.Prob.([name,'_'])=obj.make_var(dim1,obj.f);
                    else  % define up or yp parameters
                        obj.Prob.([name(1),'p_'])=obj.make_par(dim1,obj.p);
                    end
                end
                obj.Prob.rf_ = obj.make_par(obj.ny,obj.f); % define reference parameter
                er_ = obj.Prob.yf_ - obj.Prob.rf_; % error w.r.t. reference
                du_ = diff([obj.Prob.up_(:,end) obj.Prob.uf_],1,2); % u_{k+1}-u_k
                if obj.options.SolverFramework == 2 % CasADi: make user-defined constraints
                    for k = 1:length(obj.Prob.con_usr)
                       obj.Prob.Opti.subject_to(obj.Prob.con_usr{k});
                    end
                end

                % construct cost function
                obj.Prob.cost = er_(:).'*Q*er_(:)...
                               + obj.Prob.uf_(:).'*R*obj.Prob.uf_(:)...
                               + du_(:).'*dR*du_(:);
                if obj.options.SolverFramework == 2
                    obj.Prob.Opti.minimize(obj.Prob.cost);
                end
                
                % make constraints governing dynamics
                if obj.options.ExplicitPredictor
                    obj.make_con_dyn_ExplicitPredictor(); % <-- also defines G
                else
                    obj.make_con_dyn_4Optimization();
                end
                if obj.options.SolverFramework == 2
                    for k = 1:length(obj.Prob.con_dyn)
                        obj.Prob.Opti.subject_to(obj.Prob.con_dyn{k});
                    end
                end
    
                % construct optimization problem
                if solve_type.UseOptimizer
                    obj.makeOptimizer();
                    obj.solve = @obj.optimizer_solve;
                else
                    obj.solve = @obj.optimize_solve;
                end
            else
                % no user-defined constraints -> use analytical solution
                col1 = zeros(1,obj.f*obj.nu); col1(1) = 1; col1(obj.nu+1) = -1;
                row1 = zeros(1,obj.nu*obj.f); row1(1) = 1;
                Sdel = toeplitz(col1,row1);

                % cost function: J = 0.5 uf.'*H*uf + c.'*uf
                    % c = c_u0*u0 + Gu.'*Q*(Lu*up+Ly*yp-rf)
                obj.Prob.c_u0 = -Sdel.'*R(:,1:obj.nu);
                    % H = H_const + Gu.'*Q*Gu
                obj.Prob.H_const = R + Sdel.'*dR*Sdel;
                
                obj.solve = @obj.analytical_solve;
            end

            % adaptive or non-adaptive -> determines method to update data
            if options.adaptive
                obj.step_data_update = @obj.step_data_update_adaptive;
            else
                obj.step_data_update = @obj.step_data_update_non_adaptive;
            end

        end
        
        %% getter & setter functions
        function value = get.up(obj)
            value = obj.Upast(:,end-obj.p+1:end);
        end
        function set.up(obj,value)
            obj.Upast(:,end-obj.p+1:end) = value;
        end
        function value = get.yp(obj)
            value = obj.Ypast(:,end-obj.p+1:end);
        end
        function set.yp(obj,value)
            obj.Ypast(:,end-obj.p+1:end) = value;
        end

        function a = get.LHS(obj)
            % use first range 1:Nbar for compatibility with non-adaptive case
            H_u = mat2cell(obj.Upast(:,1:obj.Nbar),obj.nu,ones(1,obj.Nbar));
            H_u = H_u(hankel(1:obj.pfid,obj.pfid:obj.Nbar));
            H_y = mat2cell(obj.Ypast(:,1:obj.Nbar),obj.ny,ones(1,obj.Nbar));
            H_y = H_y(hankel(1:obj.pfid,obj.pfid:obj.Nbar));
            a = [cell2mat(H_u);cell2mat(H_y)];
            a = a/sqrt(size(a,2));
            if obj.options.use_IV
                a = a*a(1:end-obj.fid*obj.ny,:).';
            end
        end

        %% help functions

        % ============ make constraints governing dynamics ================
        function make_con_dyn_4Optimization(obj)
            % dynamics are defined by equation of the from: LHS * G = Hf
            if obj.options.SolverFramework == 1 % use YALMIP
                obj.Prob.Hf_= ...
                    [obj.make_sdp_Hankel([obj.Prob.up_ obj.Prob.uf_],obj.pfid,obj.nGcols);...
                    obj.make_sdp_Hankel([obj.Prob.yp_ obj.Prob.yf_],obj.pfid,obj.nGcols)];
            else % use CasADi
                obj.Prob.Hf_= ...
                    [obj.make_CasADi_Hankel([obj.Prob.up_ obj.Prob.uf_],obj.pfid,obj.nGcols,'u');...
                    obj.make_CasADi_Hankel([obj.Prob.yp_ obj.Prob.yf_],obj.pfid,obj.nGcols,'y')];
            end
            % define G & LHS matrix
            m1 = obj.pfid*obj.nu + obj.p*obj.ny;
            m2 = m1 + obj.fid*obj.ny;
            if obj.options.use_IV
                obj.Prob.G_   = obj.make_var(m1, obj.nGcols);
                obj.Prob.LHS_ = obj.make_par(m2, m1);
            else
                obj.Prob.G_   = obj.make_var(obj.N, obj.nGcols);
                obj.Prob.LHS_ = obj.make_par(m2,    obj.N);
            end
            
            % make constraints governing dynamics
            if obj.options.SolverFramework == 1
                obj.Prob.con_dyn = obj.Prob.LHS_*obj.Prob.G_==obj.Prob.Hf_;
            else
                obj.Prob.con_dyn = cell(1,obj.nGcols);
                for con_num = 1:obj.nGcols
                    obj.Prob.con_dyn{con_num} = obj.Prob.LHS_*obj.Prob.G_(:,con_num)==obj.Prob.Hf_(:,con_num);
                end
            end
        end
        
        %% step
        function [uf, yf_hat] = step(obj,u_k,y_k,opt)
            arguments
                obj
                u_k (:,1) double
                y_k (:,1) double
                opt.rf (:,:) double = []
            end
            % set reference if specified
            if ~isempty(opt.rf)
                obj.rf = opt.rf;
            end

            % update data
            obj.step_data_update(u_k,y_k)
            
%             if rank(obj.LHS(1:end-obj.ny,:))<min(size(obj.LHS(1:end-obj.ny,:)))
%               disp('not full rank')
%             else
%               disp('full rank')
%             end

            % solve optimization problem
            [uf, yf_hat] = obj.solve();
        end
        
        %% updating data
        function step_data_update_adaptive(obj,u_k,y_k)
            arguments
                obj
                u_k (:,1) double
                y_k (:,1) double
            end
            obj.Ypast = circshift(obj.Ypast,[0,-1]); obj.Ypast(:,end) = y_k;
            obj.Upast = circshift(obj.Upast,[0,-1]); obj.Upast(:,end) = u_k;
        end
        function step_data_update_non_adaptive(obj,u_k,y_k)
            arguments
                obj
                u_k (:,1) double
                y_k (:,1) double
            end
            obj.yp = circshift(obj.yp,[0,-1]); obj.yp(:,end) = y_k;
            obj.up = circshift(obj.up,[0,-1]); obj.up(:,end) = u_k;
        end
        
        %% optimization
        function makeOptimizer(obj)
            if obj.options.SolverFramework == 1 % use YALMIP
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
            else % use CasADi
                if obj.options.ExplicitPredictor
                    obj.Prob.Optimizer = obj.Prob.Opti.to_function('Optimizer',...
                        {obj.Prob.Lu_,obj.Prob.Ly_,obj.Prob.Gu_,obj.Prob.up_,obj.Prob.yp_,obj.Prob.rf_},... Parameters
                        {obj.Prob.uf_,obj.Prob.yf_}); % Outputs
                else
                    obj.Prob.Optimizer = obj.Prob.Opti.to_function('Optimizer',...
                        {obj.Prob.LHS_,obj.Prob.up_,obj.Prob.yp_,obj.Prob.rf_},... Parameters
                        {obj.Prob.uf_,obj.Prob.yf_});  % Outputs
                end
            end
        end
        
        % ==================== solve using 'optimizer' ====================
        function [uf, yf_hat] = optimizer_solve(obj,opt)
            % solve problem using call to yalmip Optimizer object
            arguments
                obj
                opt.rf (:,:) double = []
            end
            if ~isempty(opt.rf)
                obj.rf = opt.rf;
            end
            if obj.options.SolverFramework == 1 % using YALMIP
                try
                    if ~obj.options.ExplicitPredictor
                        [sol,errorcode] = obj.Prob.Optimizer(obj.LHS,obj.up,obj.yp,obj.rf);
                    else
                        [Lu,Ly,Gu] = obj.getPredictorMatrices();
                        [sol,errorcode] = obj.Prob.Optimizer(Lu,Ly,Gu,obj.up,obj.yp,obj.rf);
                    end
                    if errorcode~=0 && (any(isnan(sol{1})) || any(isnan(sol{2})))
                        error(yalmiperror(errorcode))
                    end
                    [uf, yf_hat] = deal(sol{:});
                catch Error
                    error(Error)
                end
            else % using CasADi
                try
                    if ~obj.options.ExplicitPredictor
                        [uf,yf_hat] = obj.Prob.Optimizer(obj.LHS,obj.up,obj.yp,obj.rf);
                    else
                        [Lu,Ly,Gu] = obj.getPredictorMatrices();
                        [uf, yf_hat] = obj.Prob.Optimizer(Lu,Ly,Gu,obj.up,obj.yp,obj.rf);
                    end
                    uf = full(uf);
                    yf_hat = full(yf_hat);
                catch Error
                    error(Error)
                end
            end
        end
        
        % ===================== solve using 'optimize' ====================
        function [uf, yf_hat] = optimize_solve(obj,opt)
            % solve problem using regular yalmip call to 'optimize'
            arguments
                obj
                opt.rf (:,:) double = []
            end
            if ~isempty(opt.rf)
                obj.rf = opt.rf;
            end
            
            if obj.options.SolverFramework == 1
                % fill in LHS matrix to enable solution with mosek
                if ~obj.options.ExplicitPredictor
                    con_dyn = replace(obj.Prob.con_dyn,obj.Prob.LHS_,obj.LHS);
                else
                    [Lu,Ly,Gu] = obj.getPredictorMatrices();
                    con_dyn = replace(obj.Prob.con_dyn,obj.Prob.Lu_,Lu);
                    con_dyn = replace(con_dyn,obj.Prob.Ly_,Ly);
                    con_dyn = replace(con_dyn,obj.Prob.Gu_,Gu);
                end
                constraints = [con_dyn;...
                               obj.Prob.con_usr;...
                               obj.Prob.rf_==obj.rf;...
                               obj.Prob.up_==obj.up;...
                               obj.Prob.yp_==obj.yp];
    
                % solve optimization problem
                diagnostics = optimize(constraints,obj.Prob.cost,obj.Prob.sdp_opts);
                if diagnostics.problem ~= 0
                    disp(yalmiperror(diagnostics.problem))
                end
                
                uf     = value(obj.Prob.uf_);
                yf_hat = value(obj.Prob.yf_);

            else % use CasADi
                if obj.options.ExplicitPredictor
                    [Lu,Ly,Gu] = obj.getPredictorMatrices();
                    obj.Prob.Opti.set_value(obj.Prob.Lu_,Lu);
                    obj.Prob.Opti.set_value(obj.Prob.Ly_,Ly);
                    obj.Prob.Opti.set_value(obj.Prob.Gu_,Gu);
                else
                    obj.Prob.Opti.set_value(obj.Prob.LHS_,obj.LHS);
                end
                obj.Prob.Opti.set_value(obj.Prob.up_,obj.up);
                obj.Prob.Opti.set_value(obj.Prob.yp_,obj.yp);
                obj.Prob.Opti.set_value(obj.Prob.rf_,obj.rf);
                sol = obj.Prob.Opti.solve();
                uf = sol.value(obj.Prob.uf_);
                yf_hat = sol.value(obj.Prob.yf_);
            end
        end

    end

    methods(Static)
        function Hankel = make_sdp_Hankel(sdp_var,dim1,dim2)
            % construct RHS for equality governing dynamics
            sdp_dim1 = size(sdp_var,1);

            % -> construct RHS for u
            sdp_cell = cell(sdp_dim1,1);
            Hankel = sdpvar(sdp_dim1*dim1,dim2,'full');
            for k_var = 1:sdp_dim1
                sdp_cell{k_var,1} = hankel(sdp_var(k_var,1:dim1),sdp_var(k_var,dim1:dim1-1+dim2));
                Hankel(k_var:sdp_dim1:end,:) = sdp_cell{k_var,1};
            end
        end

        function Hankel = make_CasADi_Hankel(CasADi_var,dim1,dim2,u_or_y)
            % construct RHS for equality governing dynamics
            var_dim1 = size(CasADi_var,1);

            % -> construct RHS block-row by block-row
            X = casadi.MX.sym('X',var_dim1,dim1+dim2-1);
            ys = cell(dim1,1);
            for i=1:dim1
                ys{i,1} = X(:,i:i+dim2-1);
            end
            Y = vertcat(ys{:});
            F = casadi.Function(['F',u_or_y],{X},{Y}); % to ensure distinct function names
            Hankel = F(CasADi_var);
        end
    end
end

