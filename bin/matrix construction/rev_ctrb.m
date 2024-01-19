function rev_ctrb_mat = rev_ctrb(s,A,B)
% creates the reversed extended observability matrix from A & B with s
% block columns

%% check inputs
% check dimensions of A and B
nx = size(A,2);
if size(B,1)~=nx
    error('Number of columns of A must be equal to the number of rows of B')
elseif size(A,1)~=nx
    error('A must be a square matrix')
end

% check validity of s - the number of block rows
if s<1 || rem(s,1)~=0
    error(['s=',num2str(s),' is invalid. Choose an integer number of block' ...
        'rows of at least 1.'])
end

%% create extended observability matrix
nu = size(B,2);
rev_ctrb_mat = zeros(nx,nu*s);
ApowerB = B;
if s > 1
    col_e = nu*(s+1);
    for col_b = nu*(s-1)+1:-nu:nu+1
        col_e = col_e-nu;
        rev_ctrb_mat(:,col_b:col_e) = ApowerB;
        ApowerB = A*ApowerB;
    end
end
rev_ctrb_mat(:,1:nu) = ApowerB;
end