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
[EUcorr_avg_all, EYcorr_avg_all] = deal(cell(1,2));
for k_e = 1:num_e
% ---------------------------- load results -------------------------------
% inputs and outputs
% u_OL = results.u_OL{k_e,1}(:,end-Nbar+1:end);
u_CL = results.u_CL(k_e,:);
% y_OL  = results.y_OL{1,1}(:,end-Nbar+1:end);
y_CL = results.y_CL(k_e,:);
[y_all, u_all] = deal(cell(1,2));
for k_c = 1:2 % loop over controllers
    u_all{k_c} = u_CL{k_c};%[u_OL u_CL{k_c}];
    y_all{k_c} = y_CL{k_c};%[y_OL y_CL{k_c}];
end
% innovation noise
% e_OL  = results.noise{k_e}(:,end-(CL_sim_steps+Nbar)+1:end-CL_sim_steps);
e_CL  = results.noise{k_e}(:,end-CL_sim_steps+1:end);
e_all = e_CL;%[e_OL e_CL];

% ------------------------ loop over controllers ---------------------------
[EUcorr_avg_run,EYcorr_avg_run,Upfid,Efid,Yp,EYcorr,EUcorr] = deal(cell(1,2));
parfor k_c = 1:2
% ------------------------ loop over time steps ---------------------------
for k1 = Nbar+1:CL_sim_steps
    k2 = k1-Nbar;
    if k2 == 1
        % get initial block-hankel matrices
        u_use = u_all{k_c}(:,k2:k1-1);
        y_use = y_all{k_c}(:,k2:k1-1);
        e_use = e_all(:,k2:k1-1);
        [Efid{k_c},Upfid{k_c},Yp{k_c}]=get_HankelMats(e_use,u_use,y_use,p,s{k_c},nu,ny);

        % calculate initial correlations
        EYcorr{k_c} = Efid{k_c}*Yp{k_c}.'/N{k_c};
        EUcorr{k_c} = Efid{k_c}*Upfid{k_c}.'/N{k_c};
        
        EYcorr_avg_run{k_c} = EYcorr{k_c};
        EUcorr_avg_run{k_c} = EUcorr{k_c};
    else
        % update block-hankel matrices
        % - noise:
        Efid{k_c}         = circshift(Efid{k_c},-1,2);
        Efid{k_c}(:,end)  = circshift(Efid{k_c}(:,end-1),-ny,1);
        Efid{k_c}(end-ny+1:end,end) = e_CL(:,k2-1);
        % - outputs:
        Yp{k_c}           = circshift(Yp{k_c},-1,2);
        Yp{k_c}(:,end)    = circshift(Yp{k_c}(:,end-1),-ny,1);
        Yp{k_c}(end-ny+1:end,end)   = y_all{k_c}(:,k1-fid{k_c}-1);
        % - inputs:
        Upfid{k_c}        = circshift(Upfid{k_c},-1,2);
        Upfid{k_c}(:,end) = circshift(Upfid{k_c}(:,end-1),-nu,1);
        Upfid{k_c}(end-nu+1:end,end)= u_CL{k_c}(:,k2-1);

        % update correlations
        if k_c == 1
            % EY
            EYcorr{k_c}(1:end-ny,1:end-ny) = EYcorr{k_c}(ny+1:end,ny+1:end);
            EYcorr{k_c}(end-ny+1:end,:)    = Efid{k_c}(end-ny+1:end,:)*Yp{k_c}.'/N{k_c};
            EYcorr{k_c}(1:end-ny,end-ny+1:end) = Efid{k_c}(1:end-ny,:)*Yp{k_c}(end-ny+1:end,:).'/N{k_c};

            % EU
            EUcorr{k_c}(1:end-ny,1:end-nu) = EUcorr{k_c}(ny+1:end,nu+1:end);
            EUcorr{k_c}(end-ny+1:end,:)    = Efid{k_c}(end-ny+1:end,:)*Upfid{k_c}.'/N{k_c};
            EUcorr{k_c}(1:end-ny,end-nu+1:end) = Efid{k_c}(1:end-ny,:)*Upfid{k_c}(end-nu+1:end,:).'/N{k_c};
        else
            % calculate correlations
            EYcorr{k_c} = Efid{k_c}*Yp{k_c}.'/N{k_c};
            EUcorr{k_c} = Efid{k_c}*Upfid{k_c}.'/N{k_c};           
        end

        % update average correlations per run
        EUcorr_avg_run{k_c} = EUcorr_avg_run{k_c}*(k2-1)/k2+EUcorr{k_c}/k2;
        EYcorr_avg_run{k_c} = EYcorr_avg_run{k_c}*(k2-1)/k2+EYcorr{k_c}/k2;
    end
    
    % plotting correlations
    %     figure(k_c)
    %     subplot(1,2,1);
    %     imagesc(EUcorr_avg_run{k_c}); set(gca, 'YDir','reverse'); colormap(gca,Cmap); caxis(gca,[-cmaxu cmaxu]); colorbar(gca);
    %     subplot(1,2,2)
    %     imagesc(EYcorr_avg_run{k_c}); set(gca, 'YDir','reverse'); colormap(gca,Cmap); caxis(gca,[-cmaxy cmaxy]); colorbar(gca);
end

if k_e == 1
    EUcorr_avg_all{k_c} = EUcorr_avg_run{k_c};
    EYcorr_avg_all{k_c} = EYcorr_avg_run{k_c};
else
    EUcorr_avg_all{k_c} = EUcorr_avg_all{k_c}*(k_e-1)/k_e+EUcorr_avg_run{k_c}/k_e;
    EYcorr_avg_all{k_c} = EYcorr_avg_all{k_c}*(k_e-1)/k_e+EYcorr_avg_run{k_c}/k_e;
end

% plotting correlations
% figure(k_c)
% subplot(1,2,1);
% imagesc(EUcorr_avg_all{k_c}); set(gca, 'YDir','reverse'); colormap(gca,Cmap); maxval1=max(abs(EUcorr_avg_all{k_c}),[],'all'); caxis(gca,[-maxval1 maxval1]);
% colorbar(gca);
% subplot(1,2,2)
% imagesc(EYcorr_avg_all{k_c}); set(gca, 'YDir','reverse'); colormap(gca,Cmap); maxval2=max(abs(EYcorr_avg_all{k_c}),[],'all'); caxis(gca,[-maxval2 maxval2]);
% colorbar(gca);

end % of loop over k_c
disp(append('Done with k_e=',num2str(k_e)))
end % of loop over k_e

%% plotting correlations
maxval1=max(cellfun(@(x) max(abs(x),[],'all'),EUcorr_avg_all),[],'all');
% maxval2=max(cellfun(@(x) max(abs(x),[],'all'),EYcorr_avg_all),[],'all');
fig5=figure('Units', 'pixels', 'pos', [80 80 1080 800],'color','white','Visible', 'on');
set(fig5,'Units','centimeters');
pos5 = get(fig5,'Position');
width5 = 8.4*1.5;
scale5 = width5/pos5(3);
set(fig5,'Position',[pos5(1:2),width5,scale5*pos5(4)])
TL = tiledlayout(fig5,1,2,'TileSpacing','compact');
ax = cell(1,2);
for k_c = 1:2
    if k_c == 1
        method_str = 'DeePC: $f_\mathrm{ID}=f$';
    else
        method_str = 'CL-DeePC: $f_\mathrm{ID}=1$';
    end
    ax{k_c}=nexttile;
    imagesc(EUcorr_avg_all{k_c}); set(gca, 'YDir','reverse'); colormap(gca,Cmap); caxis(gca,[-maxval1 maxval1]); if k_c==2; colorbar(gca);end
    title(method_str,'Interpreter','latex','FontSize',12);
    ylabel('Row index','Interpreter','latex','FontSize',12); xlabel('Column index','Interpreter','latex','FontSize',12);
    grid on; if k_c==2; yticks(1); end
    xline(ax{k_c},p*nu+0.5,'LineWidth',1.5);
%     subplot(1,2,2)
%     imagesc(EYcorr_avg_all{k_c}); set(gca, 'YDir','reverse'); colormap(gca,Cmap); caxis(gca,[-maxval2 maxval2]); colorbar(gca);
%     title({'Correlation matrix', append('$E_{i_p,',f_str,',N}Y_{i,p,N}^\top/N$')},'Interpreter','latex','FontSize',12);
%     ylabel('Row index','Interpreter','latex','FontSize',12); xlabel('Column index','Interpreter','latex','FontSize',12);
end
sgtitle(fig5,{'Average noise-input correlation matrix',append('$E_{i_p,f_\mathrm{ID},N}\left[U_{i,p,N}^\top\;\Big| U_{i_p,f_\mathrm{ID},N}^\top\right]/N$, $p=f=',num2str(p),'$, $\bar{N}=',num2str(Nbar),'$')},'Interpreter','latex','FontSize',12);
%% helper functions
function [Efid,Upfid,Yp]=get_HankelMats(eCL,uCL,yCL,p,s,nu,ny) %s{k_c},N{k_c}
    Nbar = length(eCL);
    N = Nbar-s+1;
    
    % initial hankel matrix - Ef
    H_e  = mat2cell(eCL,ny,ones(1,Nbar));
    H_e  = H_e(hankel(p+1:s,s:Nbar));
    Efid = cell2mat(H_e);
    
    % initial hankel matrix - Yp
    H_y = mat2cell(yCL,ny,ones(1,Nbar));
    H_y = H_y(hankel(1:p,p:p+N-1));
    Yp = cell2mat(H_y);
    
    % initial hankel matrix - Upf
    H_u = mat2cell(uCL,nu,ones(1,Nbar));
    H_u = H_u(hankel(1:s,s:Nbar));
    Upfid=cell2mat(H_u);
end




