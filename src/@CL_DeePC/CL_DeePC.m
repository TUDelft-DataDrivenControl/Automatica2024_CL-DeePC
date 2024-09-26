classdef CL_DeePC < Generalized_DeePC
    %CL_DEEPC Summary of this class goes here
    %   Detailed explanation goes here
    methods
        function obj = CL_DeePC(u,y,p,f,N,Q,R,dR,options,con_user,CL_opts)
            arguments
                u (:,:) double
                y (:,:) double
                p (1,1) double {mustBeInteger,mustBePositive}
                f (1,1) double {mustBeInteger,mustBePositive}
                N (1,1) double {mustBeInteger,mustBePositive}
                Q (:,:) double
                R (:,:) double
                dR (:,:) double
                options.use_IV logical   = true
                options.adaptive logical = true
                options.useAnalytic logical = true  % use analytic solution if there are no constraints
                options.ExplicitPredictor logical = true;
                options.UseOptimizer logical = true
                options.opts = []
                options.RunMakeSolver logical = true
                con_user.constr struct = struct('expr',[],'u0',[],'uf',[],'y0',[],'yf',[]);
                CL_opts.EstimateD logical = true;
            end
            fid = 1;
            options = namedargs2cell(options);
            con_user= namedargs2cell(con_user);
            obj = obj@Generalized_DeePC(u,y,p,f,fid,N,Q,R,dR,options{:}, con_user{:});
            
            % extra options for CL DeePC
            obj.options.EstimateD = CL_opts.EstimateD;
        end
        
        % ============ make constraints governing dynamics ================
        % using an explicit predictor: parameterized by 1st block-column of Gu      
        function make_con_dyn_ExplicitPredictor(obj)
            obj.Prob.Lu_ = obj.make_par(obj.f*obj.ny,obj.p*obj.nu,'Lu_');
            obj.Prob.Ly_ = obj.make_par(obj.f*obj.ny,obj.p*obj.ny,'Ly_');
            obj.Prob.Gu_ = obj.make_par(obj.f*obj.ny,obj.nu,'Gu_'); %only 1st block-column

            % make complete Gu
            if obj.options.Framework == 1 % using YALMIP
                Gu_all_ = sdpvar(obj.f*obj.ny,obj.f*obj.nu,'full');
                shape_Gu = kron(tril(toeplitz(ones(1,obj.f),1:obj.f)),ones(obj.ny,obj.nu));
                Gu_all_ = Gu_all_.*shape_Gu;
                Gu_all_(:,1:obj.nu) = obj.Prob.Gu_;
                for kc = 2:obj.f
                    r1 = obj.ny*(kc-1)+1;
                    c1 = obj.nu*(kc-1)+1;
                    c2 = obj.nu*kc;
                    Gu_all_(r1:end,c1:c2)=obj.Prob.Gu_(1:end-obj.ny*(kc-1),:);
                end
            else % using CasADi
                Gu_c1 = casadi.MX.sym('Gu1',obj.f*obj.ny,obj.nu);
                Gu = cell(1,obj.f);
                Gu{1} = Gu_c1;
                for kc = 2:obj.f
                    Gu{kc}=[zeros((kc-1)*obj.ny,obj.nu);Gu_c1(1:end-obj.ny*(kc-1),:)];
                end
                Gu_all_ = horzcat(Gu{:});

                % turn into function
                obj.Prob.GuBlkCol2Full = casadi.Function('Gu',{Gu_c1},{Gu_all_},{'Gu_BlkCol'},{'Gu_all'});
                clear Gu_c1 Gu Gu_all_;
                Gu_all_ = obj.Prob.GuBlkCol2Full(obj.Prob.Gu_);
            end
            obj.Prob.con_dyn = obj.Prob.yf_(:) == obj.Prob.Lu_*obj.Prob.up_(:) + ...
                                               obj.Prob.Ly_*obj.Prob.yp_(:) +...
                                               Gu_all_*obj.Prob.uf_(:);
            if obj.options.Framework == 2
                obj.Prob.con_dyn = {obj.Prob.con_dyn};
            end
        end

        % ================ get explicit predictor matrices ================
        function [Lu,Ly,Gu] = getPredictorMatrices(obj,varargin)
            narginchk(1,2)
            EntireGu = false; %default: only 1st block-column of Gu
            if nargin == 2
                EntireGu = varargin{1};
            end
            
            % Estimate Markov Parameters
            Up = obj.Up;
            Uf = obj.Uf;
            Yp = obj.Yp;
            Yf = obj.Yf;
            if obj.options.EstimateD % <- D also estmated implicitly
                Z = [Up;Uf;Yp];
                PredMarkov = Yf*pinv(Z); % if IV is used, same result if PE s.t. Z is full row rank
            else
                Z = [Up;Yp];
%                 if obj.options.use_IV
%                     Ziv = [Up;Uf;Yp];
%                     PredMarkov = Yf*Ziv.'*pinv(Z*Ziv.'); % Z ~= Ziv
%                 else
%                     PredMarkov = Yf*pinv(Z);
%                 end
                PredMarkov = Yf*pinv(Z);
                % add D = 0 ---> PredMarkov = [C*tKpu D C*tKpy]
                PredMarkov = [PredMarkov(:,1:obj.nu*obj.p) zeros(obj.ny,obj.nu) PredMarkov(:,obj.nu*obj.p+1:end)];
            end
            tBeta  = PredMarkov(:,1:obj.nu*(obj.p+1));      % input  predictor Markov Parameters (including D)
            tTheta = PredMarkov(:,obj.nu*(obj.p+1)+1:end);  % output predictor Markov Parameters

%             LHS_temp = obj.LHS;
%             % implicit estimation of Predictor Markov Parameters
%             At = LHS_temp(1:end-obj.ny,:);
%             Bt = LHS_temp(end-obj.ny+1:end,:);
%             tBetaTheta = Bt*pinv(At);
%             %             cond(At)
%             %             tBetaTheta = Bt/At;
%             %             tBetaTheta = lsqr(At.',Bt.').';
%             %             tBetaTheta = lsqminnorm(At.',Bt.','nowarn').';
%             %             [tBetaTheta,~] = linsolve(At.',Bt.'); tBetaTheta = tBetaTheta.';
%             tBeta  = tBetaTheta(:,1:obj.nu*(obj.p+1));      % input  predictor Markov Parameters (including D)
%             tTheta = tBetaTheta(:,obj.nu*(obj.p+1)+1:end);  % output predictor Markov Parameters
    
            % make tHf
            shape_tHf = toeplitz(obj.p+1:-1:obj.p+2-obj.f,[obj.p+1,(obj.p+2)*ones(1,obj.f-1)]);
            shape_tHf(shape_tHf<=0) = obj.p+2;
            tTheta2 = [-tTheta eye(obj.ny) zeros(obj.ny)];
            tTheta2 = mat2cell(tTheta2,obj.ny,obj.ny*ones(1,obj.p+2));
            tHf = cell2mat(tTheta2(shape_tHf));
    
            % make tLuGu & tLy
            rows_old = 1:obj.ny; % initial row numbers
            % initialize tLuGu
            n_uchan = obj.nu*(obj.p+1);
            tLuGu   = zeros(obj.ny*obj.f,n_uchan);
            tLuGu(rows_old,:) = tBeta;
            % initialize tLy
            n_ychan = obj.ny*obj.p;
            tLy   = zeros(obj.ny*obj.f,n_ychan);
            tLy(rows_old,:) = tTheta;
            for kr = 2:obj.f
                % update new row range
                rows_new = rows_old(1)+obj.ny:rows_old(end)+obj.ny;
    
                % get new rows of tLuGu
                tLuGu(rows_new,:) = circshift(tLuGu(rows_old,:),[0,obj.nu]);
                tLuGu(rows_new,1:obj.nu) = zeros(obj.ny,obj.nu);
    
                % get new rows of tLy
                tLy(rows_new,:) = circshift(tLy(rows_old,:),[0,obj.ny]);
                tLy(rows_new,1:obj.ny) = zeros(obj.ny);
    
                % update old row range
                rows_old = rows_new;
            end
            tLuGuLy = [tLuGu tLy];
    
            % tHf * LuGuLy = tLuGuLy
%             LuGuLy = fixed.forwardSubstitute(tHf.',tLuGuLy); % use forward subst. to solve
            LuGuLy = tHf\tLuGuLy;
            Lu = LuGuLy(:,1:obj.nu*obj.p);
            Gu = LuGuLy(:,obj.nu*obj.p+1:n_uchan); %only first block column of Gu
            Ly = LuGuLy(:,n_uchan+1:end);
            
            if EntireGu
                % reconstruct entire Gu from first block column
                GuCol2 = mat2cell([Gu;zeros(obj.ny,obj.nu)],obj.ny*ones(1,obj.f+1),obj.nu);
                shape_Gu = toeplitz(1:obj.f,[1 (obj.f+1)*ones(1,obj.f-1)]);
                Gu = cell2mat(GuCol2(shape_Gu));
            end
        end

        % ======================= solve analytically ======================
        function [uf, yf_hat] = analytical_solve(obj,opt)
            arguments
                obj
                opt.rf (:,:) double = []
            end
            if ~isempty(opt.rf)
                obj.rf = opt.rf;
            end
            
            [Lu,Ly,Gu] = obj.getPredictorMatrices("EntireGu",true);
             % cost function: J = 0.5 uf.'*H*uf + c.'*uf
                    % c = c_u0*u0 + Gu.'*Q*(Lu*up+Ly*yp-rf)
                    % H = H_const + Gu.'*Q*Gu
            H = obj.Prob.H_const + Gu.'*obj.Prob.Q*Gu;
            c = obj.Prob.c_u0*obj.up(:,end) + Gu.'*obj.Prob.Q*(Lu*obj.up(:)+Ly*obj.yp(:)-obj.rf(:));
            uf = -H\c;
            yf_hat = Lu*obj.up(:) + Ly*obj.yp(:) + Gu*uf;

            uf     = reshape(uf,obj.nu,obj.f);
            yf_hat = reshape(yf_hat,obj.ny,obj.f);
        end
    end
end

