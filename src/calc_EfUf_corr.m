function results = calc_EfUf_corr(results,k_e)

f = results.f;
p = results.p;
Nbar = results.Nbar;
CL_sim_steps = size(results.u_CL{1},2);
nu = size(results.u_OL{1},1);
ny = size(results.y_OL{1},1);
num_e = size(results.noise,1);
results.EfUfnorm = cell(1,2);


for k_c = 1:2
    % determine number of block-rows and columns used
    switch k_c
        case 1 % regular DeePC
            N = Nbar-(p+f-1);   % number of columns        
            fid = f;
        case 2 % CL-DeePC
            N = Nbar-p;
            fid = 1;
    end
    s = p+fid;                  % number of block rows
%     for k_e = 1:num_e
        
        % get relevant inputs and innovation noise sequences
        u_OL = results.u_OL{k_e}(:,end-Nbar+1:end);
        u_CL = results.u_CL{k_e};
        e_OL = results.noise{k_e}(:,end-(CL_sim_steps+Nbar)+1:end-CL_sim_steps);
        e_CL = results.noise{k_e}(:,end-CL_sim_steps+1:end);
        
        % initial hankel matrix - inputs
        H_u = mat2cell(u_OL,nu,ones(1,Nbar));
        H_u = H_u(hankel(p+1:s,s:Nbar));
        Ufid = cell2mat(H_u);
        u1 = Ufid(:,1);
    
        % initial hankel matrix - noise
        H_e = mat2cell(e_OL,ny,ones(1,Nbar));
        H_e = H_e(hankel(p+1:s,s:Nbar));
        Efid = cell2mat(H_e);
        e1 = Efid(:,1);
    
        % compute correlation
        EUcorr = Efid*Ufid.';
        
        % initialize norm of correlation
        EUcorr_norm    = nan(1,CL_sim_steps);
        EUcorr_norm(1) = norm(EUcorr/N,'fro')^2/fid;
        
        for k = 2:CL_sim_steps
            % update Hankel matrices
            Efid = circshift(Efid,-1,2); Efid(:,end) = circshift(Efid(:,end-1),-ny,1); Efid(end-ny+1:end,end) = e_CL(:,k-1);
            Ufid = circshift(Ufid,-1,2); Ufid(:,end) = circshift(Ufid(:,end-1),-nu,1); Ufid(end-nu+1:end,end) = u_CL(:,k-1);
            
            % new last column
            eN = Efid(:,end);
            uN = Ufid(:,end);
            
            % calculate correlation
            EUcorr = EUcorr - e1*u1.' + eN*uN.'; %Efid*Ufid.'
            EUcorr_norm(k) = norm(EUcorr/N,'fro')^2/fid;
            
            % update first column
            e1 = Efid(:,1);
            u1 = Ufid(:,1);
        end
        
        results.EfUfnorm{k_e,k_c} = EUcorr_norm;
%     end
end

end