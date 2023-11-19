function make_CasADi_Opti_solver(obj,usr_con)
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
        zero_nr = zeros(numel(lhs),1);
        
        % set new upper/lower bounds
        if strcmp(con_gleq{k},'==')
            obj.Prob.Opti.subject_to( lhs(:) == zero_nr);
        elseif strcmp(con_gleq{k},'<=') || strcmp(con_gleq{k},'=<')
            obj.Prob.Opti.subject_to( lhs(:) <= zero_nr)
        elseif strcmp(con_gleq{k},'>=') || strcmp(con_gleq{k},'=>')
            obj.Prob.Opti.subject_to( lhs(:) >= zero_nr)
        else
            error('Constraint specified incorrectly.')
        end
    end
end
%% process usr_con str2 fields
for k=1:length(str2)
    fn = str2{k};
    if isfield(usr_con,fn)
        switch k
            case 1 % u_min
                obj.Prob.Opti.subject_to( obj.Prob.uf_(:) >= repmat(usr_con.u_min,obj.f,1) );
            case 2 % u_max
                obj.Prob.Opti.subject_to( obj.Prob.uf_(:) <= repmat(usr_con.u_max,obj.f,1) );
            case 3 % y_min
                obj.Prob.Opti.subject_to( obj.Prob.yf_(:) >= repmat(usr_con.y_min,obj.f,1) );
            case 4 % y_max
                obj.Prob.Opti.subject_to( obj.Prob.yf_(:) <= repmat(usr_con.y_max,obj.f,1) );
            case 5 % du_max
                du_ = [obj.Prob.uf_(:,1)-obj.Prob.up_(:,end) obj.Prob.uf_(:,2:end)-obj.Prob.uf_(:,1:end-1)];
                du_max = repmat(usr_con.du_max,obj.f,1);
                obj.Prob.Opti.subject_to( {du_(:) >= -du_max, du_(:) <= du_max} );
            case 6 % dy_max
                dy_ = [obj.Prob.yf_(:,1)-obj.Prob.yp_(:,end) obj.Prob.yf_(:,2:end)-obj.Prob.yf_(:,1:end-1)];
                dy_max = repmat(usr_con.dy_max,obj.f,1);
                obj.Prob.Opti.subject_to( {dy_(:) >= -dy_max, dy_(:) <= dy_max} );
        end
    end
end

%% QP cost function components

% casadi symbols
H = hessian(obj.Prob.cost,obj.Prob.x_);
c = jacobian(obj.Prob.cost-1/2*obj.Prob.x_.'*H*obj.Prob.x_,obj.Prob.x_).';
c = casadi.substitute(c,obj.Prob.x_,zero_x);
obj.Prob.H = H;
obj.Prob.c = c;

%% Constraints - Dynamics
if obj.options.ExplicitPredictor
    la_new = obj.Prob.Lu_*obj.Prob.up_(:)+obj.Prob.Ly_*obj.Prob.yp_(:);
    A_new = [-obj.Prob.Gu_ speye(obj.ny*obj.f)];
else
    Hf_= [obj.make_CasADi_Hankel([obj.Prob.up_ obj.Prob.uf_],obj.pfid,obj.nGcols,'u');...
          obj.make_CasADi_Hankel([obj.Prob.yp_ obj.Prob.yf_],obj.pfid,obj.nGcols,'y')];
    x_small = [obj.Prob.uf_(:);obj.Prob.yf_(:)]; % x_ without G
    A_new = [-jacobian(Hf_(:),x_small),kron(speye(obj.nGcols),obj.Prob.LHS_)];
    la_new = casadi.substitute(Hf_(:),x_small,zeros(size(x_small)));
end
obj.Prob.Opti.subject_to( A_new*obj.Prob.x_==la_new );

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

% specify cost
obj.Prob.Opti.minimize(obj.Prob.cost);

% specify solver
obj.Prob.Opti.solver(obj.Prob.cas_opts.solver,obj.Prob.cas_opts.options);

% specify solver method - using optimizer or construction every time?
if obj.options.UseOptimizer
    if obj.options.ExplicitPredictor
        obj.Prob.Optimizer = obj.Prob.Opti.to_function('Optimizer',...
            {obj.Prob.Lu_,obj.Prob.Ly_,obj.Prob.Gu_,obj.Prob.up_,obj.Prob.yp_,obj.Prob.rf_},... Parameters
            {obj.Prob.uf_,obj.Prob.yf_}); % Outputs
    else
        obj.Prob.Optimizer = obj.Prob.Opti.to_function('Optimizer',...
            {obj.Prob.LHS_,obj.Prob.up_,obj.Prob.yp_,obj.Prob.rf_},... Parameters
            {obj.Prob.uf_,obj.Prob.yf_});  % Outputs
    end
    obj.solve = @obj.optimizer_solve;
else
    obj.solve = @obj.optimize_solve;
end