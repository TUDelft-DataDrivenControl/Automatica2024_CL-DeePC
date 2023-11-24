function make_CasADi_solver(obj,usr_con)
% Makes the solver object

str2 = {'u_min','u_max','y_min','y_max','du_max','dy_max'};
lbx = []; ubx = []; A = []; lba = []; uba = [];

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
    obj.Prob.p2Gu = casadi.Function('p2Gu',{obj.Prob.p_},{obj.Prob.Gu_});
    obj.Prob.p2Ly = casadi.Function('p2Ly',{obj.Prob.p_},{obj.Prob.Ly_});
    obj.Prob.p2Lu = casadi.Function('p2Lu',{obj.Prob.p_},{obj.Prob.Lu_});
else
    if obj.options.use_IV
        m1 = obj.pfid*obj.nu + obj.p*obj.ny;
    else
        m1 = obj.N;
    end
    m2 = m1 + obj.fid*obj.ny;
    obj.Prob.LHS_ = obj.make_par(m2,m1,'LHS');
    obj.Prob.p_ = [obj.Prob.up_(:); obj.Prob.yp_(:); obj.Prob.rf_(:); obj.Prob.LHS_(:)];
    obj.Prob.get_p = casadi.Function('get_p',... function to optain parameters
        {obj.Prob.up_,obj.Prob.yp_,obj.Prob.rf_,obj.Prob.LHS_},... parameters
        {obj.Prob.p_},... vector with parameters
        {'up','yp','rf','LHS'},{'p'}); % naming
end
zero_p = zeros(size(obj.Prob.p_));

%% optimization variables - uf, yf, G -> x
if obj.options.ExplicitPredictor
    obj.Prob.x_ = [obj.Prob.uf_(:); obj.Prob.yf_(:)];
else
    obj.Prob.G_ = obj.make_par(m1,obj.nGcols,'G');
    obj.Prob.x_ = [obj.Prob.uf_(:); obj.Prob.yf_(:); obj.Prob.G_(:)];
    
    % make functions to get results
    obj.Prob.x2G = @(x) reshape(x((obj.nu+obj.ny)*obj.f+1:(obj.nu+obj.ny)*obj.f+m1*obj.nGcols),m1,obj.nGcols);
end
obj.Prob.x2uf = @(x) reshape(x(1:obj.nu*obj.f),obj.nu,obj.f);
obj.Prob.x2yf = @(x) reshape(x(obj.nu*obj.f+1:(obj.nu+obj.ny)*obj.f),obj.ny,obj.f);
zero_x = zeros(size(obj.Prob.x_));

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

%% remove yf from optimization variable -> integrate into cost, constraints
% clear obj.Prob.yf_
% yf_past = obj.Prob.Lu_*obj.Prob.up_(:)+obj.Prob.Ly_*obj.Prob.yp_(:);
% obj.Prob.yf_ = reshape(obj.Prob.Gu_*obj.Prob.uf_(:)+yf_past,obj.ny,obj.f);
% obj.Prob.x2yf = @(x,p) reshape(...
%      obj.Prob.p2Gu(p)*x(1:obj.nu*obj.f,1)...
%     +obj.Prob.p2Lu(p)*p(1:obj.nu*obj.p)...
%     +obj.Prob.p2Ly(p)*p(obj.nu*obj.p+1:obj.p*(obj.nu+obj.ny)),...
%     obj.ny,obj.f);
% 
% % adjust matrices & bounds for lba <= A <= uba, lbx <= x <= ubx
% idx_yf  = obj.f*obj.nu+1:obj.f*(obj.nu+obj.ny); % indexes van yf in X
% mask_yf = 1:length(obj.Prob.x_);
% mask_yf = (mask_yf>obj.f*obj.nu) & (mask_yf<= obj.f*(obj.nu+obj.ny)); % true/false location of yf in X
% 
% % adjust optimization vector x
% obj.Prob.x_ = obj.Prob.x_(1:obj.nu*obj.f);
% zero_x = zeros(size(obj.Prob.x_));
% 
% % lbx & ubx
% lbx_yf = lbx(idx_yf,1);
% ubx_yf = ubx(idx_yf,1);
% lbx = lbx(~mask_yf);
% ubx = ubx(~mask_yf);
% 
% % A matrix
% Ayf = A(:,idx_yf);
% Auf = A(:,1:obj.f*obj.nu);
% Auf = Auf+Ayf*obj.Prob.Gu_;
% A = [Auf;obj.Prob.Gu_];
% 
% % lba & uba
% lba = lba-Ayf*yf_past;
% uba = uba-Ayf*yf_past;
% lba = [lba;lbx_yf-yf_past];
% uba = [uba;ubx_yf-yf_past];

er_ = obj.Prob.yf_ - obj.Prob.rf_; % error w.r.t. reference
du_ = horzcat(obj.Prob.uf_(:,1)    -obj.Prob.up_(:,end), ...
              obj.Prob.uf_(:,2:end)-obj.Prob.uf_(:,1:end-1)); % u_{k+1}-u_k
obj.Prob.cost =   er_(:).'*obj.Prob.Q *er_(:) ...
       + obj.Prob.uf_(:).'*obj.Prob.R *obj.Prob.uf_(:) ...
                + du_(:).'*obj.Prob.dR*du_(:);

%% QP cost function components
% 
% % casadi symbols
% H = hessian(obj.Prob.cost,obj.Prob.x_);
% c = jacobian(obj.Prob.cost-1/2*obj.Prob.x_.'*H*obj.Prob.x_,obj.Prob.x_).';
% c = casadi.substitute(c,obj.Prob.x_,zero_x);
% H = casadi.DM(H);
% obj.Prob.H = H;
% obj.Prob.c = c;

%% Constraints - Dynamics
if obj.options.ExplicitPredictor
    ulba_dyn = obj.Prob.Lu_*obj.Prob.up_(:)+obj.Prob.Ly_*obj.Prob.yp_(:);
    Adyn = [-obj.Prob.Gu_ speye(obj.ny*obj.f)];
else
    Hf_= [obj.make_CasADi_Hankel([obj.Prob.up_ obj.Prob.uf_],obj.pfid,obj.nGcols,'u');...
          obj.make_CasADi_Hankel([obj.Prob.yp_ obj.Prob.yf_],obj.pfid,obj.nGcols,'y')];
    x_small = [obj.Prob.uf_(:);obj.Prob.yf_(:)]; % x_ without G
    Adyn = [casadi.DM(-jacobian(Hf_(:),x_small)),kron(speye(obj.nGcols),obj.Prob.LHS_) zeros(numel(Hf_),num_S)];
    ulba_dyn = casadi.substitute(Hf_(:),x_small,zeros(size(x_small)));
end

%% make functions for regular solver
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

% regular problem
prob    = struct('f', obj.Prob.cost, 'x', obj.Prob.x_, 'g', [A;Adyn]*obj.Prob.x_,'p',obj.Prob.p_);
obj.Prob.QPsolver = casadi.qpsol('solver',obj.Prob.cas_opts.solver,prob,obj.Prob.cas_opts.options);

% make get functions
obj.Prob.get_lba = casadi.Function('get_lba',{obj.Prob.p_},{[lba;ulba_dyn]});
obj.Prob.get_uba = casadi.Function('get_uba',{obj.Prob.p_},{[uba;ulba_dyn]});
obj.Prob.get_lbx = casadi.Function('get_lbx',{obj.Prob.p_},{lbx});
obj.Prob.get_ubx = casadi.Function('get_ubx',{obj.Prob.p_},{ubx});

% initializing result structure
zero_g = zeros(size(A,1)+size(Adyn,1),1);
obj.Prob.res = struct;
obj.Prob.res.x     = zero_x;
obj.Prob.res.lam_x = zero_x;
obj.Prob.res.g     = zero_g;
obj.Prob.res.lam_g = zero_g;
obj.Prob.get_x0     = @() obj.Prob.res.x;
obj.Prob.get_lam_x0 = @() obj.Prob.res.lam_x;
obj.Prob.get_lam_a0 = @() obj.Prob.res.lam_g;

obj.Prob.p2res = @(p) obj.Prob.QPsolver('p',p,...
    'x0',obj.Prob.get_x0(),...
    'lam_x0',obj.Prob.get_lam_x0(),'lam_g0',obj.Prob.get_lam_a0(),...
    'lbg',obj.Prob.get_lba(p),'ubg',obj.Prob.get_uba(p),...
    'lbx',obj.Prob.get_lbx(p),'ubx',obj.Prob.get_ubx(p));

obj.Prob.res2ufyf = @(res) deal(full(obj.Prob.x2uf(res.x)),full(obj.Prob.x2yf(res.x)));

% specify solver method
obj.solve = @obj.optimizer_solve;

%% create backup solver for when QP is infeasible
obj.Prob.backup = struct;

% turning constraints soft (except for dynamics)
mask_ulbx = ~isinf(lbx) | ~isinf(ubx);
nxLarge = size(obj.Prob.x_,1);
P = speye(nxLarge,nxLarge); P = P(mask_ulbx,:);
na = size(A,1);
np = sum(mask_ulbx);
npa= na+np;
A = [[P;A] speye(npa,npa)];
lba = [lbx(mask_ulbx);lba];
uba = [ubx(mask_ulbx);uba];
ubx(mask_ulbx) =  Inf(np,1); ubx = [ubx; Inf(npa,1)];
lbx(mask_ulbx) = -Inf(np,1); lbx = [lbx;-Inf(npa,1)];

obj.Prob.backup.sigma_ = obj.make_var(npa,1,'sigma');
obj.Prob.backup.x_     = vertcat(obj.Prob.x_,obj.Prob.backup.sigma_);
zero_x = zeros(size(obj.Prob.backup.x_));

if obj.options.ExplicitPredictor
    Adyn = [Adyn,sparse(obj.ny*obj.f,npa)];
else
    %TODO
end

obj.Prob.backup.cost = obj.Prob.cost+1e15*obj.Prob.backup.sigma_.'*obj.Prob.backup.sigma_;

prob2    = struct('f', obj.Prob.backup.cost, 'x', obj.Prob.backup.x_, 'g', [A;Adyn]*obj.Prob.backup.x_,'p',obj.Prob.p_);
obj.Prob.backup.QPsolver = casadi.qpsol('solver',obj.Prob.cas_opts.solver,prob2,obj.Prob.cas_opts.options);

% make get functions
obj.Prob.backup.get_lba = casadi.Function('get_lba',{obj.Prob.p_},{[lba;ulba_dyn]});
obj.Prob.backup.get_uba = casadi.Function('get_uba',{obj.Prob.p_},{[uba;ulba_dyn]});
obj.Prob.backup.get_lbx = casadi.Function('get_lbx',{obj.Prob.p_},{lbx});
obj.Prob.backup.get_ubx = casadi.Function('get_ubx',{obj.Prob.p_},{ubx});

% initializing result structure
zero_g = zeros(size(A,1)+size(Adyn,1),1);
obj.Prob.backup.res = struct;
obj.Prob.backup.res.x     = zero_x;
obj.Prob.backup.res.lam_x = zero_x;
obj.Prob.backup.res.g     = zero_g;
obj.Prob.backup.res.lam_g = zero_g;
obj.Prob.backup.get_x0     = @() obj.Prob.backup.res.x;
obj.Prob.backup.get_lam_x0 = @() obj.Prob.backup.res.lam_x;
obj.Prob.backup.get_lam_a0 = @() obj.Prob.backup.res.lam_g;

obj.Prob.backup.p2res = @(p) obj.Prob.backup.QPsolver('p',p,...
    'x0',obj.Prob.backup.get_x0(),...
    'lam_x0',obj.Prob.backup.get_lam_x0(),'lam_g0',obj.Prob.backup.get_lam_a0(),...
    'lbg',obj.Prob.backup.get_lba(p),'ubg',obj.Prob.backup.get_uba(p),...
    'lbx',obj.Prob.backup.get_lbx(p),'ubx',obj.Prob.backup.get_ubx(p));
