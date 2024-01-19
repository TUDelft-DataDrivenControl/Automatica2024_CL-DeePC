function toep_mat = create_toeplitz(s,A,B,C,D)
% creates a block-toeplitz matrix with s block-rows that looks like
%   |  D                              |
%   |  CB         D                   |
%   |  CAB        CB       D          |
%   |C(A^2)B      CAB      CB  D      |
%   |   *          *       *   *  *   |
%   |C(A^s-2)B  C(A^s-3)B  *   *  CB D|

%% check inputs
% check dimensions of A, B, C, D
ny = size(C,1); % number of outputs
nu = size(B,2); % number of inputs
if ny~=size(D,1)
    error('C and D must have the same number of rows')
elseif size(A,1)~=size(A,2)
    error('A must be a square matrix')
elseif size(A,1)~=size(C,2)
    error('Number of columns of C must be equal to the number of rows of A')
elseif size(A,2)~=size(B,1)
    error('Number of columns of A must be equal to the number of rows of B')
elseif size(D,2)~=nu
    error('B and D must have an equal number of columns')
end

% check validity of s - the number of block rows
if s<1 || rem(s,1)~=0
    error(['s=',num2str(s),' is invalid. Choose an integer number of block' ...
        'rows of at least 1.'])
end

%% making block-toeplitz matrix
blocks = cell(1,s+1);
blocks{1} = D;
ApowerB = B;
if s > 2
    for row = 2:s-1
        blocks{row} = C*ApowerB;
        ApowerB = A*ApowerB;
    end
end
blocks{s} = C*ApowerB;
blocks{s+1} = zeros(ny,nu);
c = 1:s;
r = (s+1)*ones(s,1); r(1) = c(1);
toep_mat = cell2mat(blocks(toeplitz(c,r)));
end