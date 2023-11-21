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

% --------------------------------- Check flags ---------------------------
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

%% Determining solver/optimization framework
% YALMIP = 1, CasADi with Opti = 2, CasADi without Opti = 3

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

% specify how to make parameters/variables
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

%% Construct solver

if obj.options.Framework == 1
    obj.make_YALMIP_solver(usr_con,expr_flag)
elseif obj.options.Framework == 2
    obj.make_CasADi_Opti_solver(usr_con)
else
    obj.make_CasADi_solver(usr_con)
end

end