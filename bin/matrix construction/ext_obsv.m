function ext_obsv_mat = ext_obsv(s,A,C)
% creates the extended observability matrix from A & C with s block rows

%% check inputs
% check dimensions of A and C
nx = size(A,1);
if size(C,2)~=nx
    error('Number of columns of C must be equal to the number of rows of A')
elseif nx~=size(A,2)
    error('A must be a square matrix')
end

% check validity of s - the number of block rows
if s<1 || rem(s,1)~=0
    error(['s=',num2str(s),' is invalid. Choose an integer number of block' ...
        'rows of at least 1.'])
end

%% create extended observability matrix
ny = size(C,1);
ext_obsv_mat = zeros(ny*s,nx);
CApower = C;
if s > 1
    row_e = 0;
    for row_b = 1:ny:ny*(s-1)
        row_e = row_e+ny;
        ext_obsv_mat(row_b:row_e,:) = CApower;
        CApower = CApower*A;
    end
end
ext_obsv_mat(end-ny+1:end,:) = CApower;
end