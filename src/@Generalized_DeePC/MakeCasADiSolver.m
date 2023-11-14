function MakeCasADiSolver(obj,usr_con)
% Makes the solver object


%% cost function

% casadi symbols
er_ = obj.Prob.yf_ - obj.Prob.rf_; % error w.r.t. reference
du_ = [obj.Prob.uf_(:,1)-obj.Prob.up_(:,end) obj.Prob.uf_(:,2:end)-obj.Prob.uf_(:,1:end-1)]; % u_{k+1}-u_k
cost =            er_(:).'*obj.Prob.Q *er_(:) ...
       + obj.Prob.uf_(:).'*obj.Prob.R *obj.Prob.uf_(:) ...
                + du_(:).'*obj.Prob.dR*du_(:);
H = hessian(cost,obj.Prob.x_);
c = jacobian(cost-casadi.DM(1/2)*obj.Prob.x_.'*H*obj.Prob.x_,obj.Prob.x_).';
c = casadi.substitute(c,obj.Prob.x_,casadi.DM(zeros(size(obj.Prob.x_))));
H = casadi.DM(H);
obj.Prob.H = H;
obj.Prob.c = c;

%% User defined constraints
% recognize input structure fields & convert to standard used fields
fns = fields(usr_con);
for k = 1:length(fns)
    usr_fn = fns{k};
    fn1 = regexp(usr_fn,'^(d{0,1})(u|y).*(min|max)$','tokens');
    if ~isempty(fn1)
        % convert recognized field if it is different from standard
        fn1 = fn1{1};
        fn = strcat(fn1{1:2},'_',fn1{end});
        
        % get dimensions to verify
        val = usr_con.(usr_fn);
        if strcmp(fn1{2},'u')
            ndim = obj.nu;
        else
            ndim = obj.ny;
        end

        % verify dimensions
        if size(val,2)~=1 || rem(ndim,size(val,1))~=0 || size(val,1) > ndim
            error(['Incorrect input size of field ', usr_fn])
        else
            val = repmat(val,ndim/size(val,1),1);
        end

        % assign value
        usr_con.(fn) = val;
        if ~strcmp(fn,usr_fn) % do not remove field if it is the same
            usr_con = rmfield(usr_con,usr_fn);
        end
    elseif isempty(fn1)
        % discard unrecognized fields
        warning(['User-defined constraint ',usr_fn,'not recognized'])
        usr_con = rmfield(usr_con,usr_fn);
    end
end
% add missing fields
possible_fns = {'u_min','u_max','y_min','y_max','du_max','dy_max'};
for k=1:length(possible_fns)
    fn = possible_fns{k};
    fn1 = regexp(fn,'^(d{0,1})(u|y).*(min|max)$','tokens');
    if ~isfield(usr_con,fn)
        fn1 = fn1{1};
        if isempty(fn1{1}) %no du or dy => set to -inf or inf
            if strcmp(fn1{3},'min')
                Sign = -1;
            else
                Sign = 1;
            end
            if strcmp(fn1{2},'u')
                ndim = obj.nu;
            else
                ndim = obj.ny;
            end
            usr_con.(fn) = Sign*inf(ndim,1);
        else % lba <= Ax <= uba not necessarily specified
            usr_con.(fn) = [];
        end
    end
end
if obj.options.ExplicitPredictor
    lbx = [repmat(usr_con.u_min,obj.f,1);repmat(usr_con.y_min,obj.f,1)];
    ubx = [repmat(usr_con.u_max,obj.f,1);repmat(usr_con.y_max,obj.f,1)];
else
    lbx = [repmat(usr_con.u_min,obj.f,1);repmat(usr_con.y_min,obj.f,1);-inf(numel(obj.Prob.G_),1)];
    ubx = [repmat(usr_con.u_max,obj.f,1);repmat(usr_con.y_max,obj.f,1); inf(numel(obj.Prob.G_),1)];
end
uba = [];
lba = [];
A = [];
if ~isempty(usr_con.du_max)
    lba = [lba;-usr_con.du_max+obj.Prob.up_(:,end); -repmat(usr_con.du_max,obj.f-1,1)];
    uba = [uba; usr_con.du_max+obj.Prob.up_(:,end);  repmat(usr_con.du_max,obj.f-1,1)];
    du_ = [obj.Prob.uf_(:,1)-obj.Prob.up_(:,end) obj.Prob.uf_(:,2:end)-obj.Prob.uf_(:,1:end-1)];
    A   = [A;jacobian(du_(:),obj.Prob.x_)];
end
if ~isempty(usr_con.dy_max)
    lba = [lba;-usr_con.dy_max+obj.Prob.yp_(:,end); -repmat(usr_con.dy_max,obj.f-1,1)];
    uba = [uba; usr_con.dy_max+obj.Prob.yp_(:,end);  repmat(usr_con.dy_max,obj.f-1,1)];
    dy_ = [obj.Prob.yf_(:,1)-obj.Prob.yp_(:,end) obj.Prob.yf_(:,2:end)-obj.Prob.yf_(:,1:end-1)];
    A   = [A;jacobian(dy_(:),obj.Prob.x_)];
end
%% Constraints - Dynamics
if obj.options.ExplicitPredictor
    ulba = obj.Prob.Lu_*obj.Prob.up_(:)+obj.Prob.Ly_*obj.Prob.yp_(:);
    lba = [lba;ulba];
    uba = [uba;ulba];
    A = [A;[-obj.Prob.Gu_ speye(obj.ny*obj.f)]];
else
    Hf_= [obj.make_CasADi_Hankel([obj.Prob.up_ obj.Prob.uf_],obj.pfid,obj.nGcols,'u');...
          obj.make_CasADi_Hankel([obj.Prob.yp_ obj.Prob.yf_],obj.pfid,obj.nGcols,'y')];
    x_small = [obj.Prob.uf_(:);obj.Prob.yf_(:)]; % x_ without G
    A = [A;[-jacobian(Hf_(:),x_small),kron(speye(obj.nGcols),obj.Prob.LHS_)]];
    la_new = casadi.substitute(Hf_(:),x_small,casadi.DM(zeros(size(x_small))));
    lba = [lba;la_new];
    uba = [uba;la_new];
end

qp = struct();
opts = obj.Prob.cas_opts.options;
% if ~UseOptimizer
qp.h = obj.Prob.H.sparsity();
qp.a = A.sparsity();
obj.Prob.conic = casadi.conic('S',obj.Prob.cas_opts.solver,qp,opts);
% make get functions
obj.Prob.get_c   = casadi.Function('get_c',  {obj.Prob.p_},{c});
obj.Prob.get_a   = casadi.Function('get_a',  {obj.Prob.p_},{A});
obj.Prob.get_lba = casadi.Function('get_lba',{obj.Prob.p_},{lba});
obj.Prob.get_uba = casadi.Function('get_uba',{obj.Prob.p_},{uba});
obj.Prob.get_lbx = casadi.Function('get_lbx',{obj.Prob.p_},{lbx});
obj.Prob.get_ubx = casadi.Function('get_ubx',{obj.Prob.p_},{ubx});
obj.Prob.p2x = @(p) obj.Prob.conic('h',obj.Prob.H,...
                       'g',  obj.Prob.get_c(p)  ,'a',  obj.Prob.get_a(p),...
                       'lba',obj.Prob.get_lba(p),'uba',obj.Prob.get_uba(p),...
                       'lbx',obj.Prob.get_lbx(p),'ubx',obj.Prob.get_ubx(p));
if obj.options.ExplicitPredictor
    obj.Prob.Optimizer = @(Lu,Ly,Gu) obj.Prob.p2x(obj.Prob.get_p(obj.up,obj.yp,obj.rf,Lu,Ly,Gu));
else
    obj.Prob.Optimizer = @() obj.Prob.p2x(obj.Prob.get_p(obj.up,obj.yp,obj.rf,obj.LHS));
end
% else
%     qp.x = obj.Prob.x_;
%     qp.f = 0.5*obj.Prob.x_.'*H*obj.Prob.x_+c.'*obj.Prob.x_;
%     qp.g = A*obj.Prob.x_;
% %     qp.lbg = lba;
% %     qp.ubg = uba;
% %     qp.lbx = lbx;
% %     qp.ubx = ubx;
%     qp.p   = obj.Prob.p_;
%     solver = casadi.qpsol('QPsolver',obj.Prob.cas_opts.solver,qp,opts);
% end