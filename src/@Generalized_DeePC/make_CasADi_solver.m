function make_CasADi_solver(obj,usr_con)
% Makes the solver object

str2 = {'u_min','u_max','y_min','y_max','du_max','dy_max'};
lbx = []; ubx = []; A = []; lba = []; uba = [];

%% optimization variables - uf, yf, G -> x
if obj.options.ExplicitPredictor
    obj.Prob.x_ = [obj.Prob.uf_(:); obj.Prob.yf_(:)];
else
    if obj.options.use_IV
        m1 = obj.pfid*obj.nu + obj.p*obj.ny;
    else
        m1 = obj.N;
    end
    obj.Prob.G_ = obj.make_par(m1,obj.nGcols,'G');
    obj.Prob.x_ = [obj.Prob.uf_(:); obj.Prob.yf_(:); obj.Prob.G_(:)];
    
    % make functions to get results
    obj.Prob.x2G = @(x) reshape(x(end-m1*obj.nGcols+1:end),m1,obj.nGcols);
end
obj.Prob.x2uf = @(x) reshape(x(1:obj.nu*obj.f),obj.nu,obj.f);
obj.Prob.x2yf = @(x) reshape(x(obj.nu*obj.f+1:(obj.nu+obj.ny)*obj.f),obj.ny,obj.f);
zero_x = zeros(size(obj.Prob.x_));

%% parameters - up, yp, rf, Lu, Ly, Gu -> p
if obj.options.ExplicitPredictor
    obj.Prob.Lu_ = obj.make_par(obj.f*obj.ny, obj.p*obj.nu,'Lu');
    obj.Prob.Ly_ = obj.make_par(obj.f*obj.ny, obj.p*obj.ny,'Ly');
    obj.Prob.Gu_ = obj.make_par(obj.f*obj.ny, obj.f*obj.nu,'Gu');
    obj.Prob.p_ = [obj.Prob.up_(:); obj.Prob.yp_(:); obj.Prob.rf_(:); obj.Prob.Lu_(:); obj.Prob.Ly_(:); obj.Prob.Gu_(:)];
    obj.Prob.get_p = casadi.Function('get_p',... function to optain parameters
        {obj.Prob.up_,obj.Prob.yp_,obj.Prob.rf_,obj.Prob.Lu_,obj.Prob.Ly_,obj.Prob.Gu_},... parameters
        {obj.Prob.p_},... vector with parameters
        {'up','yp','rf','Lu','Ly','Gu'},{'p'}); % naming
else
    m2 = m1 + obj.fid*obj.ny;
    obj.Prob.LHS_ = obj.make_par(m2,m1,'LHS');
    obj.Prob.p_ = [obj.Prob.up_(:); obj.Prob.yp_(:); obj.Prob.rf_(:); obj.Prob.LHS_(:)];
    obj.Prob.get_p = casadi.Function('get_p',... function to optain parameters
        {obj.Prob.up_,obj.Prob.yp_,obj.Prob.rf_,obj.Prob.LHS_},... parameters
        {obj.Prob.p_},... vector with parameters
        {'up','yp','rf','LHS'},{'p'}); % naming
end
zero_p = zeros(size(obj.Prob.p_));

%% process usr_con.expr
if isfield(usr_con,'expr')
    % obtain con_LHS & con_gleq from usr_con.expr
    if ~isvector(usr_con.expr) && ismatrix(usr_con.expr)
        con_LHS  = usr_con.expr(:,1);
        con_gleq = usr_con.expr(:,2);
    elseif isvector(usr_con.expr)
        con_LHS  = usr_con.expr(1:2:end); con_LHS = con_LHS(:);
        con_gleq = usr_con.expr(2:2:end); con_gleq= con_gleq(:);
    end
    
    for k = 1:numel(con_LHS)
        lhs = con_LHS{k};

        Ax_new = casadi.DM(jacobian(lhs,obj.Prob.x_));
        % simple check for identity matrix -> use lbx, ubx
        if isdiag(sparse(Ax_new)) && all(sparse(Ax_new)==1)
            use_bx = true;
        else
            use_bx = false;
        end
        
        Ap_new = casadi.DM(jacobian(lhs,obj.Prob.p_));
        const_new = casadi.substitute(lhs,obj.Prob.x_,zero_x);
        const_new = casadi.substitute(const_new,obj.Prob.p_,zero_p);
        p_terms   = Ap_new*obj.Prob.p_;
        rhs       = -p_terms-const_new;
        num_rows  = size(const_new,1);
        
        % set new upper/lower bounds
        if strcmp(con_gleq{k},'==')
            lb_new = rhs;
            ub_new = rhs;
        elseif strcmp(con_gleq{k},'<=') || strcmp(con_gleq{k},'=<')
            lb_new = -inf(num_rows,1);
            ub_new = rhs;
        elseif strcmp(con_gleq{k},'>=') || strcmp(con_gleq{k},'=>')
            lb_new = rhs;
            ub_new = inf(num_rows,1);
        else
            error('Constraint specified incorrectly.')
        end
    
        if ~use_bx
            A   = [A;  Ax_new];
            lba = [lba;lb_new];
            uba = [uba;ub_new];
        else
            lbx = [lbx;lb_new];
            ubx = [ubx;ub_new];
        end
    end
end
%% process usr_con str2 fields
% add missing fields
for k=1:length(str2)
    fn = str2{k};
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
        else % du_max, dy_max not necessarily specified
            usr_con.(fn) = [];
        end
    end
end
if obj.options.ExplicitPredictor
    lbx = [lbx; repmat(usr_con.u_min,obj.f,1); repmat(usr_con.y_min,obj.f,1)];
    ubx = [ubx; repmat(usr_con.u_max,obj.f,1); repmat(usr_con.y_max,obj.f,1)];
else
    lbx = [lbx; repmat(usr_con.u_min,obj.f,1); repmat(usr_con.y_min,obj.f,1);-inf(numel(obj.Prob.G_),1)];
    ubx = [ubx; repmat(usr_con.u_max,obj.f,1); repmat(usr_con.y_max,obj.f,1); inf(numel(obj.Prob.G_),1)];
end
if ~isempty(usr_con.du_max)
    lba = [lba; -usr_con.du_max+obj.Prob.up_(:,end); -repmat(usr_con.du_max,obj.f-1,1)];
    uba = [uba;  usr_con.du_max+obj.Prob.up_(:,end);  repmat(usr_con.du_max,obj.f-1,1)];
    du_ = [obj.Prob.uf_(:,1)-obj.Prob.up_(:,end) obj.Prob.uf_(:,2:end)-obj.Prob.uf_(:,1:end-1)];
    A   = [A;casadi.DM(jacobian(du_(:),obj.Prob.x_))];
end
if ~isempty(usr_con.dy_max)
    lba = [lba;-usr_con.dy_max+obj.Prob.yp_(:,end); -repmat(usr_con.dy_max,obj.f-1,1)];
    uba = [uba; usr_con.dy_max+obj.Prob.yp_(:,end);  repmat(usr_con.dy_max,obj.f-1,1)];
    dy_ = [obj.Prob.yf_(:,1)-obj.Prob.yp_(:,end) obj.Prob.yf_(:,2:end)-obj.Prob.yf_(:,1:end-1)];
    A   = [A;casadi.DM(jacobian(dy_(:),obj.Prob.x_))];
end

%% QP cost function components

% casadi symbols
H = hessian(obj.Prob.cost,obj.Prob.x_);
c = jacobian(obj.Prob.cost-1/2*obj.Prob.x_.'*H*obj.Prob.x_,obj.Prob.x_).';
c = casadi.substitute(c,obj.Prob.x_,zero_x);
H = casadi.DM(H);
obj.Prob.H = H;
obj.Prob.c = c;

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
    A = [A;[casadi.DM(-jacobian(Hf_(:),x_small)),kron(speye(obj.nGcols),obj.Prob.LHS_)]];
    la_new = casadi.substitute(Hf_(:),x_small,zeros(size(x_small)));
    lba = [lba;la_new];
    uba = [uba;la_new];
end

%% make functions for solver
% solver options
opts = obj.options.opts;
if isempty(opts)
    % leave obj.Prob.cas_opts as default
elseif isstruct(opts)
    if ~isfield(opts,'solver')
        opts.solver = obj.Prob.cas_opts.solver;
    end
    obj.Prob.cas_opts = opts;
else
    error('Options passed to solver are of an unrecognized type.')
end
% selected solver-specific options
if ~isfield(obj.Prob.cas_opts,'options')
    obj.Prob.cas_opts.options = struct();
end

% Make QP solver
qp = struct();
qp.h = obj.Prob.H.sparsity();
qp.a = A.sparsity();
obj.Prob.QPsolver = casadi.conic('S',obj.Prob.cas_opts.solver,qp,obj.Prob.cas_opts.options);

% prob    = struct('f', obj.Prob.cost, 'x', obj.Prob.x_, 'g', A*obj.Prob.x_,'p',obj.Prob.p_);
% % solver  = casadi.nlpsol('solver', 'ipopt', prob,ops.opts);
% obj.Prob.QPsolver = casadi.nlpsol('solver',obj.Prob.cas_opts.solver,prob,obj.Prob.cas_opts.options);

% make get functions
obj.Prob.get_c   = casadi.Function('get_c',  {obj.Prob.p_},{c});
obj.Prob.get_a   = casadi.Function('get_a',  {obj.Prob.p_},{A});
obj.Prob.get_lba = casadi.Function('get_lba',{obj.Prob.p_},{lba});
obj.Prob.get_uba = casadi.Function('get_uba',{obj.Prob.p_},{uba});
obj.Prob.get_lbx = casadi.Function('get_lbx',{obj.Prob.p_},{lbx});
obj.Prob.get_ubx = casadi.Function('get_ubx',{obj.Prob.p_},{ubx});
% obj.Prob.get_x0  = @(obj) obj.Prob.x0;
% obj.Prob.get_lam_x0 = @(obj) obj.Prob.lam_x0;
% obj.Prob.get_lam_a0 = @(obj) obj.Prob.lam_a0;
obj.Prob.p2res = @(p) obj.Prob.QPsolver('h',obj.Prob.H,...
    'g',  obj.Prob.get_c(p)  ,'a',  obj.Prob.get_a(p),...
    'lba',obj.Prob.get_lba(p),'uba',obj.Prob.get_uba(p),...
    'lbx',obj.Prob.get_lbx(p),'ubx',obj.Prob.get_ubx(p));%,...
obj.Prob.res2ufyf = @(res) deal(obj.Prob.x2uf(res.x),obj.Prob.x2yf(res.x));

% specify solver method
obj.solve = @obj.optimizer_solve;