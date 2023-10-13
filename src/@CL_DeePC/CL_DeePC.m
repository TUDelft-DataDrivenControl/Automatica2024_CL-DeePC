classdef CL_DeePC < Generalized_DeePC
    %CL_DEEPC Summary of this class goes here
    %   Detailed explanation goes here
    methods
        function obj = CL_DeePC(u,y,p,f,N,Q,R,dR,options,con_user,solve_type)
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
                con_user.constr struct = struct('expr',[],'u0_sdp',[],'uf_sdp',[],'y0_sdp',[],'yf_sdp',[]);
                solve_type.UseOptimizer logical = true
                solve_type.sdp_opts struct = sdpsettings('solver','mosek','verbose',0);
            end
            fid = 1;
            options = namedargs2cell(options);
            con_user= namedargs2cell(con_user);
            solve_type = namedargs2cell(solve_type);
            obj = obj@Generalized_DeePC(u,y,p,f,fid,N,Q,R,dR,options{:}, con_user{:}, solve_type{:});
        end
    end
end

