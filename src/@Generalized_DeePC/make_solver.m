function make_solver(obj,usr_con)

%% parsing usr_con fields
% -> renames fields & checks dimensions & checks possible classes

% two types of constraint specifications possible using below fields
% - do not rearrange below str1 or str3!
str1 = {'uf','yf','rf','u0','y0'}; %-> fns should contain expr field too if used
str2 = {'u_min','u_max','y_min','y_max','du_max','dy_max'};
str3 = {'Opti','expr'};

% get field names
fns = fields(usr_con);

% ============ for each of the entries in str1, st2, str3: ================
% -> checks for ambiguity (multiple matches -> error)
% -> delete empty fields
% -> check class/dimensions
% -> correct field naming
for k = 1:13
    if k <= 5 % str1
        fn_target  = str1{k};   % corresponding field in str1
        fn_target1 = fn_target(1);
        fn_target2 = fn_target(end);
        reg_str = strcat('^',fn_target1,'_*',fn_target2,'$'); % field to search for
    elseif k <= 11   % str2
        fn_target  = str2{k-5}; % corresponding field in str2
        if k < 10; uy_idx = 1; else; uy_idx = 2; end
        fn_target1 = fn_target(1:uy_idx);
        fn_target2 = fn_target(end-2:end);
        reg_str = strcat('^',fn_target1,'_*',fn_target2,'$'); % field to search for
    else % str3
        fn_target = str3{k-11}; % corresponding field in str3
        reg_str = append('^',fn_target,'$'); % field to search for
    end
    
    % answer of regexpression search
    reg_ans = regexp(fns,reg_str,'tokens','ignorecase');
    
    % find fields that match
    reg_match = ~cellfun(@isempty,reg_ans); % no match -> 0, match -> 1

    % -------------------- only a single match exists ---------------------
    if sum(reg_match) == 1
        fn_usr = fns{reg_match};
        
        % -------------------- delete empty fields ------------------------
        if isempty(usr_con.(fn_usr)) && ~strcmpi(fn_usr,'Opti') % 'Opti' field can be empty
            % empty field -> delete
            usr_con = rmfield(usr_con,fn_usr);
        else
            value = usr_con.(fn_usr);

            % ------------------------ str1 & str2 ------------------------
            % => fields in str1 & str2 are of an allowed class (although possibly inconsistent)
            % => fields in str1 & str2 are of the correct dimensions
            if k <= 11 % 
                % ensure correct dimensions are used
                if strcmp(fn_target1(end),'u') % <= needed dim. from u/y/r
                    uy_dim = obj.nu;
                else
                    uy_dim = obj.ny;
                end 
                if k <= 3
                    pos_classes = {'sdpvar','casadi.MX','casadi.SX'};
                    pos_size    = [uy_dim, obj.f];
                elseif k <= 5
                    pos_classes = {'sdpvar','casadi.MX','casadi.SX'};
                    pos_size    = [uy_dim, 1];
                else
                    pos_classes = {'double'};
                    pos_size    = [uy_dim, 1];
                end
                try
                    validateattributes(value,pos_classes,{'ndims',2,'size',pos_size})
                catch Error
                    disp(append('For the field recognized as ',fn_target,' the following eror was encountered:'))
                    error(Error.message)
                end
            
            % --------------------------- Opti ----------------------------
            % => .Opti field is either empty or casadi.Opti instance
            elseif k == 12
                if ~(isempty(value) || isa(value,'casadi.Opti'))
                    error('Field recognized as Opti must be either of a casadi.Opti instance or empty.')
                end
            
            % --------------------------- expr ----------------------------
            % => expr contains only 1 'type' of constraints: SX / MX / YALMIP-type (sdpvar,constraint/lmi)
            % => if expr is cell & all constraints are YALMIP-type then they are taken out of the cell
%             elseif k == 13
%                 if isa(value,'cell')
%                     expr_is_cell = true;
%                 else % turn into cell and analyze with cellfun
%                     value = {value};
%                     expr_is_cell = false;
%                 end
%                 if isfield(usr_con,'Opti') && isempty(usr_con.Opti)
%                     if isvector(value)
%                         LHS  = value(1:2:end); LHS  = LHS(:);
%                         gleq = value(2:2:end); gleq = gleq(:);
%                     else
%                         LHS = value(:,1);
%                         gleq = value(:,2);
%                     end
%                     expr_gleq   = all(cellfun(@(x) isa(x,'char'),gleq));
%                     expr_SX     = all(cellfun(@(x) isa(x,'casadi.SX'),LHS));
%                     expr_MX     = all(cellfun(@(x) isa(x,'casadi.MX'),LHS));
%                     expr_YALMIP = all(cellfun(@(x) isa(x,'sdpvar'),LHS));
%                     
%                 else
%                     expr_SX = all(cellfun(@(x) isa(x,'casadi.SX'),value));
%                     expr_MX = all(cellfun(@(x) isa(x,'casadi.MX'),value));
%                     expr_YALMIP = all(cellfun(@(x) (isa(x,'constraint') || isa(x,'lmi')),value));
%                 end
%                 if ~( expr_SX  || expr_MX || expr_YALMIP )
%                     error("Unrecognized or inconsistent constraint types specified by constraint field 'expr'.")
%                 end
%                 if ~expr_is_cell % return to previous type
%                     value = value{1};
%                 end
%                 
%                 % turn YALMIP constraints from cell -> constraints/lmi type
%                 if expr_YALMIP && expr_is_cell
%                     value = [value{:}];
%                 end

            end
            
            % -------------------- correct field naming--------------------
            if ~strcmp(fn_target,fn_usr)
                usr_con = rmfield(usr_con,fn_usr);
            end
            usr_con.(fn_target) = value;
            fns = fields(usr_con);
        end
    
    % ------------------ multiple matches exist -> error ------------------
    elseif sum(reg_match) > 1
        error(append('Multiple constraint fields found that are similar to ',fn_target,'.'))
    end

end
% --------------- check for unrecognized fields in usr_con ----------------
str_all = [str1 ,str2 ,str3];
unrec = '';
for k = 1:length(fns)
    fn = fns{k};
    if ~any(strcmp(fn,str_all))
        unrec = append(unrec,', ',fn);
    end
end
if ~isempty(unrec)
    unrec = unrec(3:end);
    error(append('Constraint field(s) ',unrec,' have not been recognized.'))
end
% -------------------------------------------------------------------------

%% create flags for expr, str1, Opti

% ------------------------ expr handling and flag -------------------------
% possible flags: none, SX, MX, YALMIP cell, YALMIP
if isfield(usr_con,'expr')
    if isa(usr_con.expr,'cell')
        if ~isvector(usr_con.expr) && ismatrix(usr_con.expr)
            con_LHS  = usr_con.expr(:,1);
            con_gleq = usr_con.expr(:,2);
        elseif isvector(usr_con.expr)
            con_LHS  = usr_con.expr(1:2:end); con_LHS = con_LHS(:);
            con_gleq = usr_con.expr(2:2:end); con_gleq= con_gleq(:);
        else
            error('Incorrect cell array structure specified for the expr field.')
        end
        if length(con_LHS) ~= length(con_gleq)
            error('Incorrect cell array structure specified for the expr field.')
        elseif ~all(cellfun(@(x) isa(x,'char'),con_gleq))
            error('Incorrect cell array structure specified for the expr field.')
        end
        if     all(cellfun(@(x) isa(x,'casadi.SX'),con_LHS))
            expr_flag = 'SX';
        elseif all(cellfun(@(x) isa(x,'casadi.MX'),con_LHS))
            expr_flag = 'MX';
        elseif all(cellfun(@(x) isa(x,'sdpvar'),con_LHS))
            expr_flag = 'YALMIP cell';
        else
            error('Incorrect cell array structure specified for the expr field.')
        end
    elseif  isa(usr_con.expr,'constraint') || isa(usr_con.expr,'lmi')
        expr_flag = 'YALMIP';
    else
        error('Invalid constraint expression supplied.')
    end
else
    expr_flag = 'none'; % -> no 'expr' field
end

% ---------------------------------- str1 flag ----------------------------
% possible flags: none, SX, MX, YALMIP
usr_str1 = intersect(fns,str1); % -> cell array with fields that are both in usr_con & str1
if isempty(usr_str1)
    str1_flag = 'none';
else
    % get classes
    str1_class = cell(1,length(usr_str1));
    for k = 1:length(usr_str1)
        str1_class{k} = class(usr_con.(usr_str1{k}));
    end
    str1_class = unique(str1_class);
    if length(str1_class) > 1
        error('Inconsistent specification of parameters and/or variables used.')
    end
    str1_class = str1_class{1};
    if     strcmp(str1_class,'casadi.SX')
        str1_flag = 'SX';
    elseif strcmp(str1_class,'casadi.MX')
        str1_flag = 'MX';
    elseif strcmp(str1_class,'sdpvar')
        str1_flag = 'YALMIP';
    else
        error('Unrecognized class of specified parameters/variables encountered.')
    end

end

% ---------------------------------- Opti flag ----------------------------
% possible flags: none, make, present
if ~isfield(usr_con,'Opti')
    opti_flag = 'none';
elseif isempty(usr_con.Opti)
    opti_flag = 'make';
    obj.Prob.Opti = casadi.Opti('conic'); % to solve a QP problem
else
    opti_flag = 'present';
    obj.Prob.Opti = usr_con.Opti;
end

if ~contains(expr_flag,str1_flag)
    error('Inconsistent use of specified parameter/variable type and expression.');
end
if strcmp(opti_flag,'present') && ~strcmp(str1_flag,'MX')
    error('When choosing to use the caradi.Opti class, use the parameter/variable method of the created instance to specify parameters/constraints.')
end
if strcmp(opti_flag,'make') && ~strcmp(str1_flag,'none')
    error('Use and specify an instance of the casadi.Opti class to make variables and parameters.')
end
if strcmp(opti_flag,'present') && ~all([usr_con.Opti.nx usr_con.Opti.np usr_con.Opti.ng] == zeros(1,3))
    error('Specify constraints via expr field.')
end
% ============================= Further checks ============================

% % ---------------- use field in str1 <-> provide expr --------------------
% if ~isempty(usr_str1) && ~isfield(usr_con,'expr')
%     error(['Constraint fields u0, uf, y0, yf, or rf have been recognized,' ...
%           ' but a corresponding expr field for the expression has not been specified.'])
% elseif isempty(usr_str1) && isfield(usr_con,'expr')
%     error('Constraint field expr for expression specified without specifying variables and constraints used.')
% end
% 
% % ---------------------- fields in str1 of same class ---------------------
% str1_classes = cellfun(@(x) class(usr_con.(x)), usr_str1,'UniformOutput',false);
% if ~all(cellfun(@(x) strcmp(x,str1_classes{1}),str1_classes))
%     error(strjoin([{'Inconsistent types of inputs detected for constraint fields corresponding to'},usr_str1{:},{'.'}]))
% end
% 
% % ---------------------- Opti -> expr all of MX type ----------------------
% if isfield(usr_con,'Opti') && isfield(usr_con,'expr') && ( (~isempty(usr_con.Opti) && ~expr_MX) || ( isempty(usr_con.Opti) && ~(expr_MX || expr_SX) ) )
%     % is a field & ( (not empty & not all MX) or (empty & not all MX or SX ) )
%     error('For the indicated CasADi Opti framework wrong constraint types are found.\n%s', ...
%           'Specify parameters/variables using casadi.Opti.parameter/variable respectively.')
% % Note: cannot distinguish between variables/parameters made using
% % casadi.Opti.parameter/variable and casadi.MX.sym(). Opti needs former.
% elseif isfield(usr_con,'Opti') && isempty(usr_con.Opti) && isfield(usr_con,'expr') && ~expr_MX
%     error('For the indicated CasADi Opti framework wrong constraint types are found.\n%s', ...
%           'Specify parameters/variables using casadi.Opti.parameter/variable respectively.')
% elseif isfield(usr_con,'Opti') && ~isempty(usr_con.Opti) && ~all([usr_con.Opti.nx usr_con.Opti.np usr_con.Opti.ng] == zeros(1,3))
%     error('Specify constraints using expr field.')
% end

%% Determining solver/optimization framework
% YALMIP = 1, CasADi with Opti = 2, CasADi without Opti = 3
% if isfield(usr_con,'Opti')
%     if isfield(usr_con,'expr') && ( (~isempty(usr_con.Opti) && expr_MX) || ( isempty(usr_con.Opti) && (expr_MX || expr_SX) ) )
%         obj.options.Framework = 2;
%     elseif ~isfield(usr_con,'expr')
%         obj.options.Framework = 2;
%     else
%         error('Unrecognized use of Opti constraint formulation')
%     end
%     if isempty(usr_con.Opti)
%         obj.Prob.Opti = casadi.Opti('conic'); % to solve a QP problem
%     else
%         obj.Prob.Opti = usr_con.Opti;
%     end
% elseif expr_YALMIP
%     obj.options.Framework = 1;
% elseif expr_MX  || expr_SX
%     obj.options.Framework = 3;
% else
%     error('Unrecognized optimization framework for optimal control problem.')
% end

if contains(str1_flag,'YALMIP')
    obj.options.Framework = 1;
    make_parvar_method = 1;
elseif strcmp(str1_flag,'MX') && strcmp(opti_flag,'present')
    obj.options.Framework = 2; % note: cannot distinguish between Opti.variable/parameter and MX.sym
    make_parvar_method = 2;
elseif strcmp(str1_flag,'MX')
    obj.options.Framework = 3;
    make_parvar_method = 3;
elseif strcmp(str1_flag,'SX')
    obj.options.Framework = 3;
    make_parvar_method = 4;
elseif strcmp(str1_flag,'none') && ~strcmp(opti_flag,'none')
    obj.options.Framework = 2;
    make_parvar_method = 2;
elseif strcmp(str1_flag,'none') &&  strcmp(opti_flag,'none') % use casadi w/ SX symbolics
    obj.options.Framework = 3;
    make_parvar_method = 4;
else
    error('Unrecognized use of constraint structure.')
end

% determine how to make parameters/variables
switch make_parvar_method
    case 1
        obj.make_var = @(dim1,dim2,name) sdpvar(dim1,dim2,'full');
        obj.make_par = @(dim1,dim2,name) sdpvar(dim1,dim2,'full');
    case 2
        obj.make_var = @(dim1,dim2,name) obj.Prob.Opti.variable(dim1,dim2);
        obj.make_par = @(dim1,dim2,name) obj.Prob.Opti.parameter(dim1,dim2);
    case 3
        obj.make_var = @(dim1,dim2,name) casadi.MX.sym(name,dim1,dim2);
        obj.make_par = @(dim1,dim2,name) casadi.MX.sym(name,dim1,dim2);
    case 4
        obj.make_var = @(dim1,dim2,name) casadi.SX.sym(name,dim1,dim2);
        obj.make_par = @(dim1,dim2,name) casadi.SX.sym(name,dim1,dim2);
end

%% make variables/parameters up, yp, rf, uf, yf

if isfield(usr_con,'u0')
    obj.Prob.up_ = [obj.make_par(obj.nu,obj.p-1,'up_endmin1'),usr_con.u0];
else
    obj.Prob.up_ = obj.make_par(obj.nu,obj.p,'up');
end
if isfield(usr_con,'y0')
    obj.Prob.yp_ = [obj.make_par(obj.ny,obj.p-1,'yp_endmin1'),usr_con.y0];
else
    obj.Prob.yp_ = obj.make_par(obj.nu,obj.p,'yp');
end
if isfield(usr_con,'rf')
    obj.Prob.rf_ = usr_con.rf;
else
    obj.Prob.rf_ = obj.make_par(obj.ny,obj.f,'rf');
end
if isfield(usr_con,'uf')
    obj.Prob.uf_ = usr_con.uf;
else
    obj.Prob.uf_ = obj.make_var(obj.nu,obj.f,'uf');
end
if isfield(usr_con,'yf')
    obj.Prob.yf_ = usr_con.yf;
else
    obj.Prob.yf_ = obj.make_var(obj.ny,obj.f,'yf');
end

% %% optimization variables - uf, yf, G -> x
% if obj.options.ExplicitPredictor
%     obj.Prob.x_ = [obj.Prob.uf_(:); obj.Prob.yf_(:)];
% else
%     if obj.options.use_IV
%         m1 = obj.pfid*obj.nu + obj.p*obj.ny;
%     else
%         m1 = obj.N;
%     end
%     obj.Prob.G_ = obj.make_par(m1,obj.nGcols,'G');
%     obj.Prob.x_ = [obj.Prob.uf_(:); obj.Prob.yf_(:); obj.Prob.G_(:)];
%     
%     % make functions to get results
%     obj.Prob.x2G = @(x) reshape(x(end-m1*obj.nGcols+1:end),m1,obj.nGcols);
% end
% obj.Prob.x2uf = @(x) reshape(x(1:obj.nu*obj.f),obj.nu,obj.f);
% obj.Prob.x2yf = @(x) reshape(x(obj.nu*obj.f+1:(obj.nu+obj.ny)*obj.f),obj.ny,obj.f);
% zero_x = zeros(size(obj.Prob.x_));
% 
% %% parameters - up, yp, rf, Lu, Ly, Gu -> p
% if obj.options.ExplicitPredictor
%     obj.Prob.Lu_ = obj.make_par(obj.f*obj.ny, obj.p*obj.nu,'Lu');
%     obj.Prob.Ly_ = obj.make_par(obj.f*obj.ny, obj.p*obj.ny,'Ly');
%     obj.Prob.Gu_ = obj.make_par(obj.f*obj.ny, obj.f*obj.nu,'Gu');
%     obj.Prob.p_ = [obj.Prob.up_(:); obj.Prob.yp_(:); obj.Prob.rf_(:); obj.Prob.Lu_(:); obj.Prob.Ly_(:); obj.Prob.Gu_(:)];
%     if obj.options.Framework > 1
%         obj.Prob.get_p = casadi.Function('get_p',... function to optain parameters
%                     {obj.Prob.up_,obj.Prob.yp_,obj.Prob.rf_,obj.Prob.Lu_,obj.Prob.Ly_,obj.Prob.Gu_},... parameters
%                     {obj.Prob.p_},... vector with parameters
%                     {'up','yp','rf','Lu','Ly','Gu'},{'p'}); % naming
%     end
% else
%     m2 = m1 + obj.fid*obj.ny;
%     obj.Prob.LHS_ = obj.make_par(m2,m1,'LHS');
%     obj.Prob.p_ = [obj.Prob.up_(:); obj.Prob.yp_(:); obj.Prob.rf_(:); obj.Prob.LHS_(:)];
%     if obj.options.Framework > 1
%         obj.Prob.get_p = casadi.Function('get_p',... function to optain parameters
%                         {obj.Prob.up_,obj.Prob.yp_,obj.Prob.rf_,obj.Prob.LHS_},... parameters
%                         {obj.Prob.p_},... vector with parameters
%                         {'up','yp','rf','LHS'},{'p'}); % naming
%     end
% end
% zero_p = zeros(size(obj.Prob.p_));

%% Cost
er_ = obj.Prob.yf_ - obj.Prob.rf_; % error w.r.t. reference
du_ = horzcat(obj.Prob.uf_(:,1)    -obj.Prob.up_(:,end), ...
              obj.Prob.uf_(:,2:end)-obj.Prob.uf_(:,1:end-1)); % u_{k+1}-u_k
obj.Prob.cost =   er_(:).'*obj.Prob.Q *er_(:) ...
       + obj.Prob.uf_(:).'*obj.Prob.R *obj.Prob.uf_(:) ...
                + du_(:).'*obj.Prob.dR*du_(:);

%% Construct solver

if obj.options.Framework == 1
    obj.make_YALMIP_solver(usr_con,expr_flag)
elseif obj.options.Framework == 2
    obj.make_CasADi_Opti_solver(usr_con)
else
    obj.make_CasADi_solver(usr_con)
end
% %% constraints specified as in str1
% 
% % initialize matrices
% uba = []; lba = []; A = []; lbx = []; ubx = [];
% 
% % ---------------------- form YALMIP constraints --------------------------
% % -> if they are still inside cell array (expr_flag = 'YALMIP cell')
% if contains(expr_flag,'YALMIP')
%     if contains(expr_flag,'cell')
%         usr_con.expr = [];
%         for k = 1:numel(con_LHS)
%             if strcmp(con_gleq{k},'==')
%                 usr_con.expr = [usr_con.expr;con_LHS{k}==0];
%             elseif strcmp(con_gleq{k},'<=') || strcmp(con_gleq{k},'=<')
%                 usr_con.expr = [usr_con.expr;con_LHS{k}<=0];
%             elseif strcmp(con_gleq{k},'>=') || strcmp(con_gleq{k},'=>')
%                 usr_con.expr = [usr_con.expr;con_LHS{k}>=0];
%             else
%                 error('Constraint specified incorrectly.')
%             end
%         end
%     end
% 
% % ---------------------- form CasADi constraints --------------------------
% elseif strcmp(expr_flag,'MX') || strcmp(expr_flag,'SX')
%     for k = 1:numel(con_LHS)
%         Ax_jac =  casadi.DM(jacobian(con_LHS{k},obj.Prob.x_));
%         Ap_jac =  casadi.DM(jacobian(con_LHS{k},obj.Prob.p_));
%         p_terms= Ap_jac*obj.Prob.p_;
%         Ac     =  casadi.substitute(con_LHS{k}, obj.Prob.x_, zero_x);
%         Ac     =  casadi.DM(casadi.substitute(Ac, obj.Prob.p_, zero_p));
%         rhs    = -p_terms-Ac;
%         num_rows = size(con_LHS{k},1);
%         A = [A;Ax_jac];
%         if strcmp(con_gleq{k},'==')
%             uba = [uba;rhs];
%             lba = [lba;rhs];
%         elseif strcmp(con_gleq{k},'<=') || strcmp(con_gleq{k},'=<')
%             uba = [uba;rhs];
%             lba = [lba;-inf(num_rows,1)];
%         elseif strcmp(con_gleq{k},'>=') || strcmp(con_gleq{k},'=>')
%             uba = [uba;inf(num_rows,1)];
%             lba = [lba;rhs];
%         else
%             error('Constraint specified incorrectly.')
%         end
%     end
% end
% %% constraints specified as in str2
% 
% % add missing fields with -inf/inf
% for k=1:length(str2)
%     fn = str2{k};
%     fn1 = regexp(fn,'^(d{0,1})(u|y).*(min|max)$','tokens');
%     if ~isfield(usr_con,fn)
%         fn1 = fn1{1};
%         if isempty(fn1{1}) %no du or dy => set to -inf or inf
%             if strcmp(fn1{3},'min')
%                 Sign = -1;
%             else
%                 Sign = 1;
%             end
%             if strcmp(fn1{2},'u')
%                 ndim = obj.nu;
%             else
%                 ndim = obj.ny;
%             end
%             usr_con.(fn) = Sign*inf(ndim,1);
%         else % lba <= Ax <= uba not necessarily specified
%             usr_con.(fn) = [];
%         end
%     end
% end
% 
% if obj.options.ExplicitPredictor
%     lbx = [lbx;repmat(usr_con.u_min,obj.f,1);repmat(usr_con.y_min,obj.f,1)];
%     ubx = [ubx;repmat(usr_con.u_max,obj.f,1);repmat(usr_con.y_max,obj.f,1)];
% else
%     lbx = [lbx;repmat(usr_con.u_min,obj.f,1);repmat(usr_con.y_min,obj.f,1);-inf(numel(obj.Prob.G_),1)];
%     ubx = [ubx;repmat(usr_con.u_max,obj.f,1);repmat(usr_con.y_max,obj.f,1); inf(numel(obj.Prob.G_),1)];
% end
% 
% % construct lba, uba, A
% if ~isempty(usr_con.du_max)
%     lba = [lba;-usr_con.du_max+obj.Prob.up_(:,end); -repmat(usr_con.du_max,obj.f-1,1)];
%     uba = [uba; usr_con.du_max+obj.Prob.up_(:,end);  repmat(usr_con.du_max,obj.f-1,1)];
%     du_ = [obj.Prob.uf_(:,1)-obj.Prob.up_(:,end) obj.Prob.uf_(:,2:end)-obj.Prob.uf_(:,1:end-1)];
%     A   = [A;jacobian(du_(:),obj.Prob.x_)];
% end
% if ~isempty(usr_con.dy_max)
%     lba = [lba;-usr_con.dy_max+obj.Prob.yp_(:,end); -repmat(usr_con.dy_max,obj.f-1,1)];
%     uba = [uba; usr_con.dy_max+obj.Prob.yp_(:,end);  repmat(usr_con.dy_max,obj.f-1,1)];
%     dy_ = [obj.Prob.yf_(:,1)-obj.Prob.yp_(:,end) obj.Prob.yf_(:,2:end)-obj.Prob.yf_(:,1:end-1)];
%     A   = [A;jacobian(dy_(:),obj.Prob.x_)];
% end
% 
% %% Constraints - Dynamics
% if obj.options.ExplicitPredictor
%     ulba = obj.Prob.Lu_*obj.Prob.up_(:)+obj.Prob.Ly_*obj.Prob.yp_(:);
%     lba = [lba;ulba];
%     uba = [uba;ulba];
%     A = [A;[-obj.Prob.Gu_ speye(obj.ny*obj.f)]];
% else
%     Hf_= [obj.make_CasADi_Hankel([obj.Prob.up_ obj.Prob.uf_],obj.pfid,obj.nGcols,'u');...
%           obj.make_CasADi_Hankel([obj.Prob.yp_ obj.Prob.yf_],obj.pfid,obj.nGcols,'y')];
%     x_small = [obj.Prob.uf_(:);obj.Prob.yf_(:)]; % x_ without G
%     A = [A;[-jacobian(Hf_(:),x_small),kron(speye(obj.nGcols),obj.Prob.LHS_)]];
%     la_new = casadi.substitute(Hf_(:),x_small,casadi.DM(zeros(size(x_small))));
%     lba = [lba;la_new];
%     uba = [uba;la_new];
% end
% 
% %% Cost function
% 
% H = hessian(cost,obj.Prob.x_);
% c = jacobian(cost-casadi.DM(1/2)*obj.Prob.x_.'*H*obj.Prob.x_,obj.Prob.x_).';
% c = casadi.substitute(c,obj.Prob.x_,zero_x);
% H = casadi.DM(H);
% obj.Prob.H = H;
% obj.Prob.c = c;
% 
% %%
% qp = struct();
% opts = obj.Prob.cas_opts.options;
% % if ~UseOptimizer
% qp.h = obj.Prob.H.sparsity();
% qp.a = A.sparsity();
% obj.Prob.conic = casadi.conic('S',obj.Prob.cas_opts.solver,qp,opts);
% % make get functions
% obj.Prob.get_c   = casadi.Function('get_c',  {obj.Prob.p_},{c});
% obj.Prob.get_a   = casadi.Function('get_a',  {obj.Prob.p_},{A});
% obj.Prob.get_lba = casadi.Function('get_lba',{obj.Prob.p_},{lba});
% obj.Prob.get_uba = casadi.Function('get_uba',{obj.Prob.p_},{uba});
% obj.Prob.get_lbx = casadi.Function('get_lbx',{obj.Prob.p_},{lbx});
% obj.Prob.get_ubx = casadi.Function('get_ubx',{obj.Prob.p_},{ubx});
% obj.Prob.p2x = @(p) obj.Prob.conic('h',obj.Prob.H,...
%                        'g',  obj.Prob.get_c(p)  ,'a',  obj.Prob.get_a(p),...
%                        'lba',obj.Prob.get_lba(p),'uba',obj.Prob.get_uba(p),...
%                        'lbx',obj.Prob.get_lbx(p),'ubx',obj.Prob.get_ubx(p));
% if obj.options.ExplicitPredictor
%     obj.Prob.Optimizer = @(Lu,Ly,Gu) obj.Prob.p2x(obj.Prob.get_p(obj.up,obj.yp,obj.rf,Lu,Ly,Gu));
% else
%     obj.Prob.Optimizer = @() obj.Prob.p2x(obj.Prob.get_p(obj.up,obj.yp,obj.rf,obj.LHS));
% end

end