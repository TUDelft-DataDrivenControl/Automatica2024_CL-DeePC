close all
clc;
% get parameters from results structure
f = results.f;
p = results.p;
nu = size(results.u_OL{1},1);
ny = size(results.y_OL{1},1);
num_e = size(results.noise,1);
Nbar = results.Nbar;
CL_sim_steps = size(results.u_CL{1},2);
ref   = results.ref(1:CL_sim_steps+f-1);
% determine dimensions of block-Hankel matrices used
N = cell(1,2);  % # of columns
s = N;          % # of block rows of u & y
N{1} = Nbar-(p+f-1); % regular DeePC
s{1} = p+f;
fid{1} = f;
N{2} = Nbar-p;       % CL-DeePC
s{2} = p+1;
fid{2} = 1;

%% create colormap
num_colors = 1001;
c1 = [220,50,32]/256;
c2 = ones(1,3);
c3 = [0,90,181]/256;
cmaxu = 2.0e-2;
cmaxy = 0.9e0;
Cmap = [interp1([0,0.5,1],[c1(1),c2(1),c3(1)],linspace(0,1,num_colors));...
        interp1([0,0.5,1],[c1(2),c2(2),c3(2)],linspace(0,1,num_colors));...
        interp1([0,0.5,1],[c1(3),c2(3),c3(3)],linspace(0,1,num_colors))].';

%% loop over k_e
[UfZp_corr_avg_all, UYcorr_avg_all] = deal(cell(1,2));
for k_e = 1%1:num_e
% ---------------------------- load results -------------------------------
% inputs and outputs
u_OL = results.u_OL{k_e,1}(:,end-Nbar+1:end);
u_CL = results.u_CL(k_e,:);
y_OL  = results.y_OL{1,1}(:,end-Nbar+1:end);
y_CL = results.y_CL(k_e,:);
[y_all, u_all] = deal(cell(1,3));
for k_c = 1:2 % loop over controllers
    u_all{k_c} = [u_OL u_CL{k_c}];
    y_all{k_c} = [y_OL y_CL{k_c}];
end

% ------------------------ loop over controllers ---------------------------
[Upfid,Yp] = deal(cell(1,2));
for k_c = 1:2
% ------------------------ loop over time steps ---------------------------
[Smin1,Smax1,cond_UfUf,cond_UfUf_Zp,Smin2,Smax2,cond_ZpZp,cond_ZpZp_Uf] = deal(nan(1,CL_sim_steps));
for k1 = Nbar+1:Nbar+CL_sim_steps
    k2 = k1-Nbar;
    if k2 == 1
        % get initial block-hankel matrices
        u_use = u_all{k_c}(:,k2:k1-1);
        y_use = y_all{k_c}(:,k2:k1-1);
        [Upfid{k_c},Yp{k_c}]=get_HankelMats(Nbar,u_use,y_use,p,s{k_c},nu,ny);
        
        % make Up, Ufid
        Up{k_c}   = Upfid{k_c}(1:p*nu,:);
        Ufid{k_c} = Upfid{k_c}(p*nu+1:end,:);

        % calculate initial correlations
        UfUp_corr{k_c} = Ufid{k_c}*Up{k_c}.'/N{k_c};
        UfYp_corr{k_c} = Ufid{k_c}*Yp{k_c}.'/N{k_c};
        UfUf_corr{k_c} = Ufid{k_c}*Ufid{k_c}.'/N{k_c};
        YpYp_corr{k_c} = Yp{k_c}*Yp{k_c}.'/N{k_c};
        UpUp_corr{k_c} = Up{k_c}*Up{k_c}.'/N{k_c};
        UpYp_corr{k_c} = Up{k_c}*Yp{k_c}.'/N{k_c};
    else
        % update block-hankel matrices
        % - outputs:
        Yp{k_c}           = circshift(Yp{k_c},-1,2);
        Yp{k_c}(:,end)    = circshift(Yp{k_c}(:,end-1),-ny,1);
        Yp{k_c}(end-ny+1:end,end)   = y_all{k_c}(:,k1-fid{k_c}-1);
        % - inputs:
        Upfid{k_c}        = circshift(Upfid{k_c},-1,2);
        Upfid{k_c}(:,end) = circshift(Upfid{k_c}(:,end-1),-nu,1);
        Upfid{k_c}(end-nu+1:end,end)= u_CL{k_c}(:,k2-1);
        
        % make Up, Ufid
        Up{k_c}   = Upfid{k_c}(1:p*nu,:);
        Ufid{k_c} = Upfid{k_c}(p*nu+1:end,:);

        % update correlations
        if k_c == 1
            UfUp_corr{k_c} = update_correlation(Ufid{k_c},Up{k_c},  UfUp_corr{k_c},nu,nu,N{k_c}); % UfUp
            UfYp_corr{k_c} = update_correlation(Ufid{k_c},Yp{k_c},  UfYp_corr{k_c},nu,ny,N{k_c}); % UfYp
            UfUf_corr{k_c} = update_correlation(Ufid{k_c},Ufid{k_c},UfUf_corr{k_c},nu,nu,N{k_c}); % UfUf
            YpYp_corr{k_c} = update_correlation(Yp{k_c},  Yp{k_c},  YpYp_corr{k_c},ny,ny,N{k_c}); % YpYp
            UpUp_corr{k_c} = update_correlation(Up{k_c},  Up{k_c},  UpUp_corr{k_c},nu,nu,N{k_c}); % UpUp
            UpYp_corr{k_c} = update_correlation(Up{k_c},  Yp{k_c},  UpYp_corr{k_c},nu,ny,N{k_c}); % UpYp
        else
            % calculate correlations
            UfUp_corr{k_c} = Ufid{k_c}*Up{k_c}.'/N{k_c};
            UfYp_corr{k_c} = Ufid{k_c}*Yp{k_c}.'/N{k_c};
            UfUf_corr{k_c} = Ufid{k_c}*Ufid{k_c}.'/N{k_c};
            YpYp_corr{k_c} = Yp{k_c}*Yp{k_c}.'/N{k_c};
            UpUp_corr{k_c} = Up{k_c}*Up{k_c}.'/N{k_c};
            UpYp_corr{k_c} = Up{k_c}*Yp{k_c}.'/N{k_c};
        end
        
%         % update average correlations per run
%         EUcorr_avg_run{k_c} = EUcorr_avg_run{k_c}*(k2-1)/k2+EUcorr{k_c}/k2;
%         EYcorr_avg_run{k_c} = EYcorr_avg_run{k_c}*(k2-1)/k2+EYcorr{k_c}/k2;
    end
    Zp = [Yp{k_c};Up{k_c}];
    ZpZp_corr = [YpYp_corr{k_c} UpYp_corr{k_c}.';UpYp_corr{k_c} UpUp_corr{k_c}];%Zp*Zp.'/N{k_c};
    UfZp_corr = [UfYp_corr{k_c} UfUp_corr{k_c}];

    % calculate Cholesky factorizations
    Chol_Uf = chol(UfUf_corr{k_c}).'; %Ufid{k_c}/sqrt(N{k_c});
    Chol_Up = chol(UpUp_corr{k_c}).';
    Chol_Yp = chol(YpYp_corr{k_c}).';
    Chol_Zp = chol(ZpZp_corr).'; %Zp/sqrt(N{k_c});
    
    normalized_UfZp_corr = Chol_Uf\UfZp_corr/(Chol_Zp.');
    normalized_UfUp_corr = Chol_Uf\UfUp_corr{k_c}/(Chol_Up.');
    normalized_UfYp_corr = Chol_Uf\UfYp_corr{k_c}/(Chol_Yp.');
    
    % plotting correlations
    figure(k_c)
    imagesc(normalized_UfZp_corr);
    set(gca, 'YDir','reverse'); colormap(gca,Cmap); colorbar(gca); caxis(gca,[-0.3 0.3]);
    
    % effect on Gu estimate
    norm_RHS1 = normalized_UfZp_corr*normalized_UfZp_corr.';
    Sall1 = eig(norm_RHS1);
    Smax1(k2) = max(Sall1,[],'all')^2;
    Smin1(k2) = min(Sall1,[],'all')^2;
    cond_UfUf(k2) = cond(UfUf_corr{k_c});
    UfUf_Zp = UfUf_corr{k_c}-UfZp_corr/ZpZp_corr*UfZp_corr.';%Chol_Uf*(eye(nu*fid{k_c})-norm_RHS1)*Chol_Uf.';
    cond_UfUf_Zp(k2) = cond(UfUf_Zp);
    % effect on Obsv*Kp estimate
    norm_RHS2 = normalized_UfZp_corr.'*normalized_UfZp_corr;
    Sall2 = eig(norm_RHS2);
    Smax2(k2) = max(Sall2,[],'all')^2;
    Smin2(k2) = min(Sall2,[],'all')^2;
    cond_ZpZp(k2) = cond(ZpZp_corr);
    ZpZp_Uf = ZpZp_corr-UfZp_corr.'/UfUf_corr{k_c}*UfZp_corr;%Chol_Uf*(eye(nu*fid{k_c})-norm_RHS1)*Chol_Uf.';
    cond_ZpZp_Uf(k2) = cond(ZpZp_Uf);
    if rem(k2,100) == 0
        figure(2+k_c)
        if k_c==1
            title_str = 'DeePC';
        else
            title_str = 'CL-DeePC';
        end
        % effect on Gu
        ax1=subplot(2,1,1);
        plot(ax1,cond_UfUf.*(1-Smin1)./(1-Smax1)); hold on;
        plot(ax1,cond_UfUf_Zp); hold off;
        ylabel('Effect on Gu'); ax1.YScale = 'log';

        % effect on Obsv*Kp
        ax2=subplot(2,1,2);
        plot(ax2,cond_ZpZp.*(1-Smin2)./(1-Smax2)); hold on;
        plot(ax2,cond_ZpZp_Uf); hold off;
        ylabel('Effect on Obsv*Kp'); ax2.YScale = 'log';
        sgtitle(title_str);
    end
end

end % of loop over k_c
disp(append('Done with k_e=',num2str(k_e)))
end % of loop over k_e

%% plotting correlations
maxval1=max(cellfun(@(x) max(abs(x),[],'all'),EUcorr_avg_all),[],'all');
% maxval2=max(cellfun(@(x) max(abs(x),[],'all'),EYcorr_avg_all),[],'all');
fig5=figure('Units', 'pixels', 'pos', [80 80 1080 480],'color','white','Visible', 'on');
TL = tiledlayout(fig5,1,2,'TileSpacing','tight');
ax = cell(1,2);
for k_c = 1:2
    if k_c == 1
        method_str = 'DeePC: $f_\mathrm{ID}=f$';
    else
        method_str = 'CL-DeePC: $f_\mathrm{ID}=1$';
    end
    ax{k_c}=nexttile;
    imagesc(EUcorr_avg_all{k_c}); set(gca, 'YDir','reverse'); colormap(gca,Cmap); caxis(gca,[-maxval1 maxval1]); colorbar(gca);
    title(method_str,'Interpreter','latex','FontSize',12);
    ylabel('Row index','Interpreter','latex','FontSize',12); xlabel('Column index','Interpreter','latex','FontSize',12);
    grid on;
    xline(ax{k_c},p*nu+0.5,'LineWidth',1.5);
%     subplot(1,2,2)
%     imagesc(EYcorr_avg_all{k_c}); set(gca, 'YDir','reverse'); colormap(gca,Cmap); caxis(gca,[-maxval2 maxval2]); colorbar(gca);
%     title({'Correlation matrix', append('$E_{i_p,',f_str,',N}Y_{i,p,N}^\top/N$')},'Interpreter','latex','FontSize',12);
%     ylabel('Row index','Interpreter','latex','FontSize',12); xlabel('Column index','Interpreter','latex','FontSize',12);
end
sgtitle(fig5,{'Average noise-input correlation matrix',append('$E_{i_p,f_\mathrm{ID},N}\left[U_{i,p,N}^\top\;\Big| U_{i_p,f_\mathrm{ID},N}^\top\right]/N$, $p=f=',num2str(p),'$, $\bar{N}=',num2str(Nbar),'$')},'Interpreter','latex','FontSize',14);
%% helper functions
function [Upfid,Yp]=get_HankelMats(Nbar,uCL,yCL,p,s,nu,ny) %s{k_c},N{k_c}
    N = Nbar-s+1;
    
    % initial hankel matrix - Yp
    H_y = mat2cell(yCL,ny,ones(1,Nbar));
    H_y = H_y(hankel(1:p,p:p+N-1));
    Yp = cell2mat(H_y);
    
    % initial hankel matrix - Upf
    H_u = mat2cell(uCL,nu,ones(1,Nbar));
    H_u = H_u(hankel(1:s,s:Nbar));
    Upfid=cell2mat(H_u);
end

function corr_new = update_correlation(H1,H2,corr_old,n1,n2,N)
    corr_new = nan(size(corr_old));
    corr_new(1:end-n1,1:end-n2) = corr_old(n1+1:end,n2+1:end);
    corr_new(end-n1+1:end,:)    = H1(end-n1+1:end,:)*H2.'/N;
    corr_new(1:end-n1,end-n2+1:end) = H1(1:end-n1,:)*H2(end-n2+1:end,:).'/N;
end