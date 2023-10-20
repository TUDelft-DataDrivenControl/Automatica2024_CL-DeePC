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
                      'con_usr',[],'con_dyn',[],'sdp_opts',[])
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
                options.ExplicitPredictor logical = false % for CL DeePC
                con_user.constr struct = struct('expr',[],'u0_sdp',[],'uf_sdp',[],'y0_sdp',[],'yf_sdp',[]);
                solve_type.UseOptimizer logical = true
                solve_type.sdp_opts struct = sdpsettings('solver','mosek','verbose',0);
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
            if ~isempty(con_user.constr.expr)
                % user-defined constraints involved -> use optimization

                % define optimization variables
                obj.Prob.yf_ = sdpvar(obj.ny,obj.f,   'full');
                obj.Prob.uf_ = sdpvar(obj.nu,obj.f,   'full');
                obj.Prob.rf_ = sdpvar(obj.ny,obj.f,   'full');
                obj.Prob.yp_ = sdpvar(obj.ny,obj.p,   'full');
                obj.Prob.up_ = sdpvar(obj.nu,obj.p,   'full');
                    % tracking error & u increments
                er_ = obj.Prob.yf_ - obj.Prob.rf_; % error w.r.t. reference
                du_ = diff([obj.Prob.up_(:,end) obj.Prob.uf_],1,2); % u_{k+1}-u_k

                % construct cost function
                obj.Prob.cost = er_(:).'*Q*er_(:)...
                               + obj.Prob.uf_(:).'*R*obj.Prob.uf_(:)...
                               + du_(:).'*dR*du_(:);
    
                % parse user-defined constraints, replace user's uf & yf
                S = namedargs2cell(con_user.constr);
                obj.Prob.con_usr = parse_user_con(obj,S{:});
                
                % make constraints governing dynamics
                if obj.options.ExplicitPredictor
                    obj.make_con_dyn_ExplicitPredictor(); % <-- also defines G
                else
                    obj.make_con_dyn_4Optimization();
                end
    
                % construct optimization problem
                obj.Prob.sdp_opts = sdpsettings(solve_type.sdp_opts);
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
            obj.Prob.Hf_= ...
                [obj.make_sdp_Hankel([obj.Prob.up_ obj.Prob.uf_],obj.pfid,obj.nGcols);...
                obj.make_sdp_Hankel([obj.Prob.yp_ obj.Prob.yf_],obj.pfid,obj.nGcols)];
            
            % define G & LHS matrix
            m1 = obj.pfid*obj.nu + obj.p*obj.ny;
            m2 = m1 + obj.fid*obj.ny;
            if obj.options.use_IV
                obj.Prob.G_   = sdpvar(m1, obj.nGcols, 'full');
                obj.Prob.LHS_ = sdpvar(m2, m1,         'full');
            else
                obj.Prob.G_   = sdpvar(obj.N, obj.nGcols, 'full');
                obj.Prob.LHS_ = sdpvar(m2,    obj.N,      'full');
            end

            % make constraints governing dynamics
            obj.Prob.con_dyn = obj.Prob.LHS_*obj.Prob.G_==obj.Prob.Hf_;
        end
        
        % ================ parse user-defined constraints =================
        function usr_con = parse_user_con(obj,usr_con)
                % Parses constraint structure. Returns user defined
                % constraints with variables defined by user switched for
                % those used in this class.
                % Field meanings:
                %  -> expr:   expressions of constraints
                %  -> uf_sdp: sdp variable used to represent uf in constraints
                %  -> u0_sdp: sdp variable used to represent last input
                %  -> yf_sdp: sdp variable used to represent yf in constraints
                %  -> y0_sdp: sdp variable used to represent last output
                arguments
                    obj
                    usr_con.uf_sdp = []
                    usr_con.yf_sdp = []
                    usr_con.u0_sdp = []
                    usr_con.y0_sdp = []
                    usr_con.expr   = []
                end
                if ~isempty(usr_con.expr)
                    struct1 = struct(... user_var, {dims, sdp variable used}
                        'uf_sdp',{[obj.nu, obj.f],obj.Prob.uf_},...
                        'yf_sdp',{[obj.ny, obj.f],obj.Prob.yf_},...
                        'u0_sdp',{[obj.nu, 1    ],obj.Prob.up_(:,end)},...
                        'y0_sdp',{[obj.ny, 1    ],obj.Prob.yp_(:,end)});
                    fns1 = fieldnames(struct1);
                    for kk = 1:length(fns1)
                        fn1 = fns1{kk};
                        user_var = usr_con.(fn1);
                        if ~isempty(user_var)
                            dims = struct1(1).(fn1);
                            check_sdp_dims(user_var,dims(1),dims(2));
                            usr_con.expr = replace(usr_con.expr,user_var,struct1(2).(fn1));
                        end
                    end
                end
                usr_con = usr_con.expr;

                function check_sdp_dims(sdpvariable,dim1,dim2)
                    % checks sdpvar size & whether it is 'full'
                    if ~all(size(sdpvariable) == [dim1,dim2])
                        [s1,s2] = size(sdpvariable);
                        str1 = ['[', num2str(s1) ,' x ', num2str(s2),']'];
                        str2 = ['[', num2str(dim1) ,' x ', num2str(dim2),']'];
                        error(['Size of sdp variable: ',str1, ' must be of size: ',str2]);
                    elseif numel(getvariables(sdpvariable)) ~= dim1*dim2
                        error('Sdp variable is not full')
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
                error(yalmiperror(errorcode))
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
    end
end

