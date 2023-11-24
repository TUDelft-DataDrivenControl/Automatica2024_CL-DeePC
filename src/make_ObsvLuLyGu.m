function [Obsv_f,Lu,Ly,Gu] = make_ObsvLuLyGu(plant,f,p)
[A,B,C,D] = ssdata(plant);
nx = size(A,1);
ny = size(C,1);
nu = size(B,2)-ny;

K = B(:,end-ny+1:end);
B = B(:,1:nu);
D = D(:,1:nu);

Atilde = A-K*C;
Btilde = B-K*D;

%% Extended observability matrix
Obsv_f = zeros(ny*f,nx);
Obsv_f(1:ny,:) = C;
for k = 2:f
    Obsv_f((k-1)*ny+1:k*ny,:) = Obsv_f((k-2)*ny+1:(k-1)*ny,:)*A;
end

%% Ktilde_pu
Ktilde_pu = zeros(nx,p*nu);
a = Btilde;
Ktilde_pu(:,end-nu+1:end) = Btilde;
for p_k = 2:p
    a = Atilde*a;
    Ktilde_pu(:,end-p_k*nu+1:end-(p_k-1)*nu) = a;
end

%% Ktilde_py
Ktilde_py = zeros(nx,p*ny);
a = K;
Ktilde_py(:,end-ny+1:end) = K;
for p_k = 2:p
    a = Atilde*a;
    Ktilde_py(:,end-p_k*ny+1:end-(p_k-1)*ny) = a;
end

%% Lu & Ly
Lu = Obsv_f*Ktilde_pu;
Ly = Obsv_f*Ktilde_py;

%% Tuf - toeplitz matrix describing influence of inputs
Tuf_struct = toeplitz(2:f+1,[2 ones(1,f-1)]);
Tuf_cells = cell(1,f+1);
Tuf_cells{1} = zeros(ny,nu);
Tuf_cells{2} = D;
Tuf_cells(3:end) = mat2cell(Obsv_f(1:end-ny,:),ny*ones(1,f-1),nx);
Tuf_cells(3:end) = cellfun(@(x) x*B,Tuf_cells(3:end),'UniformOutput',false);
Gu = cell2mat(Tuf_cells(Tuf_struct));
end