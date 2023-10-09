classdef CL_DeePC < handle
    %CL_DEEPC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        p
        f
        N
        Nbar
        nu
        ny

        % data
        Up_Yp_U1
        Y1
        up
        yp

        % optimization variables & functions
        cost
        yf
        uf
        rf
        G
    end
    
    methods
        function obj = CL_DeePC(u,y,p,f,N,Qk,Rk)
            %CL_DEEPC Construct an instance of this class
            %   Detailed explanation goes here
            obj.nu = min(size(u)); if size(u,2) < size(u,1); u=u.'; end
            obj.ny = min(size(y)); if size(y,2) < size(y,1); y=y.'; end

            % N >= n + (p+1)*nu + p*ny
            % using p as lower bound for n
            obj.p = p;
            obj.f = f;
            Nmin = ...obj.p + 
                (obj.p+1)*obj.nu + obj.p*obj.ny;
            if N < Nmin
                obj.N = Nmin;
                disp('Larger number of columns necessary for persistency of excitation')
            else
                obj.N = N;
            end
            obj.Nbar = obj.N+obj.p;

            obj.Up_Yp_U1 = [data2hankel(obj.p, u(:,1:obj.p+obj.N-1),  N);... inputs: s,data,N, use oldest [false,default] / oldest [true]?
                            data2hankel(obj.p, y(:,1:obj.p+obj.N-1),  N);
                            u(:,obj.p+1:obj.Nbar)];
%             if rank(obj.Up_Yp_U1)<min(size(obj.Up_Yp_U1))
%                 disp('not full rank')
%             else
%                 disp('full rank')
%             end
            obj.Y1       =  y(:,obj.p+1:obj.Nbar);
            obj.up       =  u(:,end-p+1:end);
            obj.yp       =  y(:,end-p+1:end);
            
            % make optimization problem
            construct_cost(obj,Qk,Rk);
        end

        function construct_cost(obj,Qk,Rk)
            obj.yf = sdpvar(obj.ny,obj.f,'full');
            obj.rf = sdpvar(obj.ny,obj.f,'full');
            obj.uf = sdpvar(obj.nu,obj.f,'full');
            
            obj.cost = (obj.yf(:)-obj.rf(:)).'*kron(eye(obj.f),Qk)*(obj.yf(:)-obj.rf(:)) + obj.uf(:).'*kron(eye(obj.f),Rk)*obj.uf(:);
        end

        function [uf, yf_hat, G] = optimize_solve(obj,rf)
%             obj.G  = sdpvar(size(obj.Up_Yp_U1,1), obj.f,'full');
            obj.G = sdpvar(obj.N, obj.f,'full');
            up_ = sdpvar(obj.nu,obj.p,'full');
            yp_ = sdpvar(obj.ny,obj.p,'full');
            
            % construct RHS for u
            u_all = [up_ obj.uf];
            u_cell = cell(obj.nu,1);
            u_RHS = sdpvar(obj.nu*(obj.p+1),obj.f,'full');
            for k_u = 1:obj.nu
                u_cell{k_u,1} = hankel(u_all(k_u,1:obj.p+1),u_all(k_u,obj.p+1:obj.p+obj.f));
                u_RHS(k_u:obj.nu:end,:) = u_cell{k_u,1};
            end

            % construct RHS for y
            y_all = [yp_ obj.yf];
            y_cell = cell(obj.ny,1);
            y_RHS = sdpvar(obj.ny*(obj.p+1),obj.f,'full');
            for k_y = 1:obj.ny
                y_cell{k_y,1} = hankel(y_all(k_y,1:obj.p+1),y_all(k_y,obj.p+1:obj.p+obj.f));
                y_RHS(k_y:obj.ny:end,:) = y_cell{k_y,1};
            end
            constraints = [obj.Up_Yp_U1;obj.Y1]*obj.G==[u_RHS;y_RHS];
%             constraints = [obj.Up_Yp_U1;obj.Y1]*obj.Up_Yp_U1.'*obj.G==[u_RHS;y_RHS];
            constraints = [constraints;...
                           obj.rf==rf;...
                           up_==obj.up;...
                           yp_==obj.yp];
%             options = sdpsettings('solver','mosek');
            diagnostics = optimize(constraints,obj.cost);%,options);
            if diagnostics.problem == 0
             disp('Solver thinks it is feasible')
            elseif diagnostics.problem == 1
             disp('Solver thinks it is infeasible')
            else
             disp('Something else happened')
            end
            
            uf     = value(obj.uf);
            yf_hat = value(obj.yf);
            G      = value(obj.G);
        end

        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

