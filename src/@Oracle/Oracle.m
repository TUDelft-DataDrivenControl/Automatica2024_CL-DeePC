classdef Oracle < handle
    %ORACLEDEEPC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Obsv_f
        Gu
        up
        yp
        rf
        nx
        nu
        ny
        p
        f
        QP = struct;
        res2uf
        res2yf
    end
    
    methods
        function obj = Oracle(Obsv_f,Gu,f,p,nx,nu,ny,Q,R,dR,con)
            % makes an Oracle solver with assumed perfect knowledge of the 
            % matrices Obsv_f, Gu, and the initial state

            [obj.Obsv_f,obj.Gu,obj.p,obj.f,obj.nx,obj.nu,obj.ny] = deal(Obsv_f,Gu,p,f,nx,nu,ny);

            %% optimization variables
            x0_ = casadi.SX.sym('x0',nx,1);
            u0_ = casadi.SX.sym('u0',nu,1);
            rf_ = casadi.SX.sym('rf',ny,f);
            uf_ = casadi.SX.sym('uf',nu,f);
            yf_ = casadi.SX.sym('yf',ny,f);
            p_ = [x0_;u0_;rf_(:)];
            x_ = [uf_(:);yf_(:)];

            %% cost
            er_f = yf_(:)-rf_(:);
            duf_ = [u0_-uf_(:,1), uf_(:,1:end-1) - uf_(:,2:end)];
            cost = er_f.'*Q*er_f + uf_(:).'*R*uf_(:) + duf_(:).'*dR*duf_(:);

            %% constraints

            % dynamics
            qdyn  = yf_(:) - Obsv_f*x0_ - Gu*uf_(:);
            qdynb = zeros(ny*f,1);

            % user defined constraints
            % du max
            q   = [qdyn;   duf_(:)];
            qlb = [qdynb; -repmat(con.du_max,f,1)];
            qub = [qdynb;  repmat(con.du_max,f,1)];

            % umax & umin
            xlb = [repmat(con.u_min,f,1);repmat(con.y_min,f,1)];
            xub = [repmat(con.u_max,f,1);repmat(con.y_max,f,1)];

            %% Make solver
            prob   = struct('f', cost, 'x', x_, 'g', q,'p',p_);
            Solver = 'ipopt';
            ipopt_opts = struct('print_time',0,'ipopt',struct('print_level',0,'nlp_scaling_method','none','warm_start_init_point','yes','hessian_constant','yes','mehrotra_algorithm','yes','max_iter',20));
            obj.QP.solver1 = casadi.nlpsol('Oracle',Solver,prob,ipopt_opts);

            % initialize solver guesses
            obj.QP.res.x      = zeros(f*(nu+ny),1);
            obj.QP.res.lam_x = zeros(f*(nu+ny),1);
            obj.QP.res.lam_g = zeros(size(qlb));

            % make functions to get result
            obj.QP.get_x0       = @() obj.QP.res.x(1:f*(nu+ny),1);      % indices specified to accomodate use after backup solver
            obj.QP.get_lam_x0   = @() obj.QP.res.lam_x(1:f*(nu+ny),1);
            obj.QP.get_lam_g0   = @() obj.QP.res.lam_g(1:size(qlb,1),1);
            obj.QP.get_res1 = @(x0_s,u0_s,ref_s) obj.QP.solver1('p',[x0_s;u0_s;ref_s(:)],...
                'x0',obj.QP.get_x0(),'lam_x0',obj.QP.get_lam_x0(),'lam_g0',obj.QP.get_lam_g0(),...
                'lbg',qlb,'ubg',qub,'lbx',xlb,'ubx',xub);
            obj.QP.res2uf   = @(res) reshape(full(res.x(1:nu*f)),nu,f);
            obj.QP.res2yf   = @(res) reshape(full(res.x(nu*f+1:(nu+ny)*f)),ny,f);

            %% Make backup solver <- solves problem with softened user-defined constraints
            
            % create slack variables
            n_sigma = (2*nu+ny)*f;
            sigma = casadi.SX.sym('sigma',n_sigma,1);
            x2_ = vertcat(x_,sigma);
            xlb2 = [xlb+sigma(1:(nu+ny)*f,1);-inf(n_sigma,1)];
            xub2 = [xub+sigma(1:(nu+ny)*f,1); inf(n_sigma,1)];
            qlb2 = [qdynb; -repmat(con.du_max,f,1)+sigma((nu+ny)*f+1:end,1)];
            qub2 = [qdynb;  repmat(con.du_max,f,1)+sigma((nu+ny)*f+1:end,1)];

            % new cost
            cost2 = cost + 1e15*(sigma.'*sigma);
            
            % create new solver
            prob2   = struct('f', cost2, 'x', x2_, 'g', q,'p',p_);
            obj.QP.solver2 = casadi.nlpsol('Oracle',Solver,prob2,ipopt_opts);
            obj.QP.get_res2 = @(x0_s,u0_s,ref_s) obj.QP.solver1('p',[x0_s;u0_s;ref_s(:)],...
                'lbg',qlb2,'ubg',qub2,'lbx',xlb2,'ubx',xub2);
        end
        
        %% Other methods
        function [uf,yf_hat,varargout] = solve(obj,x_k1,u_k0,opt)
            arguments
                obj
                x_k1 (:,1) double
                u_k0 (:,1) double
                opt.rf (:,:) double = []
            end
            if ~isempty(opt.rf)
                obj.rf = opt.rf;
            end
            % -> initial guess: uf = 0, yf = Lu*up+Ly*yp
            yf0 = obj.Obsv_f*x_k1;
            obj.QP.res.x  = [zeros(obj.nu*obj.f,1);yf0];

            try
                % try to solve regular problem
                obj.QP.stat = 0;  % <- indicates unsuccessful solve(s)
                obj.QP.res = obj.QP.get_res1(x_k1,u_k0,obj.rf);
                obj.QP.stat = 1;  % indicating solver used
            catch
                try
                    % try to solve original problem with softened constraints
                    % -> solve relaxed problem
                    obj.QP.stat = 2;
                    obj.QP.res = obj.QP.get_res2(x_k1,u_k0,obj.rf);
                    obj.QP.stat = 3;
                catch Error
                    disp(['No feasible solution found. Solution status:', num2str(obj.QP.stat) ])
                    error(Error.message);
                end
            end
            uf     = obj.QP.res2uf(obj.QP.res);
            yf_hat = obj.QP.res2yf(obj.QP.res);
            varargout{1} = obj.QP.stat; % solution status
        end
    end
end

