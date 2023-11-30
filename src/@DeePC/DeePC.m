classdef DeePC < Generalized_DeePC
    %CL_DEEPC Summary of this class goes here
    %   Detailed explanation goes here
    methods
        function obj = DeePC(u,y,p,f,N,Q,R,dR,options,con_user)
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
                con_user.constr struct = struct('expr',[],'u0',[],'uf',[],'y0',[],'yf',[]);
            end
            fid = f;
            options = namedargs2cell(options);
            con_user= namedargs2cell(con_user);
            obj = obj@Generalized_DeePC(u,y,p,f,fid,N,Q,R,dR,options{:}, con_user{:});
        end
        
        % ============ make constraints governing dynamics ================
        % using an explicit predictor: parameterized by Gu
        function make_con_dyn_ExplicitPredictor(obj)
            obj.Prob.Lu_ = obj.make_par(obj.f*obj.ny,obj.p*obj.nu);
            obj.Prob.Ly_ = obj.make_par(obj.f*obj.ny,obj.p*obj.ny);
            obj.Prob.Gu_ = obj.make_par(obj.f*obj.ny,obj.f*obj.nu);
            obj.Prob.con_dyn = obj.Prob.yf_(:) == obj.Prob.Lu_*obj.Prob.up_(:) + ...
                                               obj.Prob.Ly_*obj.Prob.yp_(:) + ...
                                               obj.Prob.Gu_*obj.Prob.uf_(:);
            if obj.options.Framework == 2
                obj.Prob.con_dyn = {obj.Prob.con_dyn};
            end
        end

        % ================ get explicit predictor matrices ================
        function [Lu,Ly,Gu] = getPredictorMatrices(obj,varargin)
                        
            % Estimate Markov Parameters
            Up = obj.Up;
            Uf = obj.Uf;
            Yp = obj.Yp;
            Yf = obj.Yf;
            Z = [Up;Uf;Yp];
            if obj.options.use_IV
                LuGuLy = Yf*Z.'*pinv(Z*Z.');
            else
                LuGuLy = Yf*pinv(Z);
            end
%             LHS_temp = obj.LHS;
% %             implicit estimation of Predictor Markov Parameters
%             At = LHS_temp(1:end-obj.ny*obj.f,:);
%             Bt = LHS_temp(end-obj.ny*obj.f+1:end,:);
%                         LuGuLy = Bt*pinv(At);
            %             cond(Z)
            %             tBetaTheta = Bt/At;
            %             tBetaTheta = lsqr(At.',Bt.').';
            %             tBetaTheta = lsqminnorm(At.',Bt.','nowarn').';
            %             [tBetaTheta,~] = linsolve(At.',Bt.'); tBetaTheta = tBetaTheta.';
            Lu = LuGuLy(:,1:obj.nu*obj.p);
            Gu = LuGuLy(:,obj.nu*obj.p+1:obj.nu*(obj.p+obj.f));
            Ly = LuGuLy(:,obj.nu*(obj.p+obj.f)+1:end);
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
            
            [Lu,Ly,Gu] = obj.getPredictorMatrices();
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

