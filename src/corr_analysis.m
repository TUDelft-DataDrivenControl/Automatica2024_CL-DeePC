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

% get data from results structure
u_OL = results.u_OL{k_e,1}(:,end-Nbar+1:end);
u_CL = results.u_CL(k_e,:);
y_OL  = results.y_OL{1,1}(:,end-Nbar+1:end);
y_CL = results.y_CL(k_e,:);
e_OL  = results.noise{k_e}(:,end-(CL_sim_steps+Nbar)+1:end-CL_sim_steps);
e_CL  = results.noise{k_e}(:,end-CL_sim_steps+1:end);
e_all = [e_OL e_CL];
u_all = cell(1,3);
y_all = u_all;
for k = 1:3
    u_all{k} = [u_OL u_CL{k}];
    y_all{k} = [y_OL y_CL{k}];
end
ref   = results.ref(1:CL_sim_steps+f-1);

results.EfUfnorm = cell(1,2);
%% Basic figure
screensize = get( 0, 'ScreenSize' );
fig1 = figure('units','pixels','position',[0 0 screensize(3:4)*0.9]);%[0 0 1440 1080]);
sgtitle(append('$\sigma^2(e_k)=',num2str(results.Re),'$, $\sigma^2(d^\mathrm{u}_k)=',num2str(results.Rdu),...
    '$, $\bar{N}=',num2str(results.Nbar),'$, $p=',num2str(results.p),'$, $f=',num2str(results.f),...
    '$, $Q=',num2str(results.Q),'$, $R=',num2str(results.R),'$, $R_\Delta=',num2str(results.dR),'$'),'Interpreter','latex','FontSize',12);
ax1=subplot(3,2,1:2);
plot(Nbar+1:Nbar+CL_sim_steps,y_CL{1}); hold on
plot(Nbar+1:Nbar+CL_sim_steps,y_CL{2});
plot(Nbar+1:Nbar+CL_sim_steps,y_CL{3});
plot(Nbar+1:Nbar+CL_sim_steps+f-1,ref);
plot(1:Nbar,y_OL,HandleVisibility='off');
xline(Nbar+0.5,HandleVisibility='off');
xline(2*Nbar+0.5,HandleVisibility='off');
legend({'DeePC','CL-DeePC','oracle','reference'},'Interpreter','latex','FontSize',12);
ylabel('$y$','Interpreter','latex','FontSize',12);
xlabel('time step','Interpreter','latex','FontSize',12);
grid on;
y_lims = ylim(ax1); y_min = y_lims(1); y_max = y_lims(2);

%% Plotting rest
% create movie frame-holder
F = repmat(getframe(fig1),1,CL_sim_steps);
vidfile = VideoWriter('testmovie3.mp4','MPEG-4');
open(vidfile);

% create colormap
num_colors = 1001;
c1 = [220,50,32]/256;
c2 = ones(1,3);
c3 = [0,90,181]/256;
cmaxu = 0.7e0;
cmaxy = 1.2e1;
Cmap = [interp1([0,0.5,1],[c1(1),c2(1),c3(1)],linspace(0,1,num_colors));...
        interp1([0,0.5,1],[c1(2),c2(2),c3(2)],linspace(0,1,num_colors));...
        interp1([0,0.5,1],[c1(3),c2(3),c3(3)],linspace(0,1,num_colors))].';

% determine dimensions of block-Hankel matrices used
N = cell(1,2);  % # of columns
s = N;          % # of block rows of u & y
N{1} = Nbar-(p+f-1); % regular DeePC
s{1} = p+f;
N{2} = Nbar-p;       % CL-DeePC
s{2} = p+1;

% subplots
ax1u = subplot(3,2,3); % regular DeePC - Ef*Upfid.'
title(ax1u,'Regular DeePC - $E_{f_\mathrm{ID}}U_{pf_\mathrm{ID}}^\top$','Interpreter','latex','FontSize',10); hold on;
xlim([0.5 s{1}*nu+0.5]); box off;
ax2u = subplot(3,2,4); % CL-DeePC      - Ef*Upfid.'
title(ax2u,'CL-DeePC - $E_{f_\mathrm{ID}}U_{pf_\mathrm{ID}}^\top$','Interpreter','latex','FontSize',10); hold on;
xlim([0.5 s{2}*nu+0.5]); box off;
ax1y = subplot(3,2,5); % regular DeePC - Ef*Yp.'
title(ax1y,'Regular DeePC - $E_{f_\mathrm{ID}}Y_p^\top$','Interpreter','latex','FontSize',10); hold on;
xlim([0.5 p*ny+0.5]); box off;
ax2y = subplot(3,2,6); % CL-DeePC      - Ef*Yp.'
title(ax2y,'CL-DeePC - $E_{f_\mathrm{ID}}Y_p^\top$','Interpreter','latex','FontSize',10); hold on;
xlim([0.5 p*ny+0.5]); box off;

fillcol = [128 128 128]/255; % fill color: grey
update_interval = 5;
[e1,u1,y1,eN,uN,yN,Upfid,Efid,Yp,EYcorr,EUcorr] = deal(cell(1,2));
for k = Nbar+1:update_interval:Nbar+CL_sim_steps
    k2 = k-Nbar;
    e_use = e_all(:,k-Nbar:k-1);
    for k_c = 1:2
        u_use = u_all{k_c}(:,k-Nbar:k-1);
        y_use = y_all{k_c}(:,k-Nbar:k-1);
        [Efid{k_c},Yp{k_c},Upfid{k_c}]=get_HankelMats(e_use,u_use,y_use,N{k_c},p,s{k_c},nu,ny);
        EYcorr{k_c} = Efid{k_c}*Yp{k_c}.';
        EUcorr{k_c} = Efid{k_c}*Upfid{k_c}.';
    %{
    if k == 1
        % initial hankel matrix - Ef
        H_e = mat2cell(e_OL,ny,ones(1,Nbar));
        H_e = H_e(hankel(p+1:s{k_c},s{k_c}:Nbar));
        Efid{k_c} = cell2mat(H_e);
        e1{k_c} = Efid{k_c}(:,1);

        % initial hankel matrix - Yp
        H_y = mat2cell(y_OL,nu,ones(1,Nbar));
        H_y = H_y(hankel(1:p,p:p+N{k_c}-1));
        Yp{k_c}  = cell2mat(H_y);
        y1{k_c} = Yp{k_c}(:,1);

        % initial hankel matrix - Upf
        H_u = mat2cell(u_OL,nu,ones(1,Nbar));
        H_u = H_u(hankel(1:s{k_c},s{k_c}:Nbar));
        Upfid{k_c}=cell2mat(H_u);
        u1{k_c} = Upfid{k_c}(:,1);
        
        % compute initial correlations based on OL data
        EYcorr{k_c} = Efid{k_c}*Yp{k_c}.';
        EUcorr{k_c} = Efid{k_c}*Upfid{k_c}.';
    else
        % update block-hankel matrices
        % - noise:
        Efid{k_c}         = circshift(Efid{k_c},-1,2);
        Efid{k_c}(:,end)  = circshift(Efid{k_c}(:,end-1),-ny,1);
        Efid{k_c}(end-ny+1:end,end) = e_CL(:,k-1);
        % - outputs:
        Yp{k_c}           = circshift(Yp{k_c},-1,2);
        Yp{k_c}(:,end)    = circshift(Yp{k_c}(:,end-1),-ny,1);
        Yp{k_c}(end-ny+1:end,end)   = y_CL{k_c}(:,k-1);
        % - inputs:
        Upfid{k_c}        = circshift(Upfid{k_c},-1,2);
        Upfid{k_c}(:,end) = circshift(Upfid{k_c}(:,end-1),-nu,1);
        Upfid{k_c}(end-nu+1:end,end)= u_CL{k_c}(:,k-1);

        % update last column
        eN{k_c} = Efid{k_c}(:,end);
        yN{k_c} = Yp{k_c}(:,end);
        uN{k_c} = Upfid{k_c}(:,end);

        % update correlations
        EYcorr{k_c} = EYcorr{k_c}-e1{k_c}*y1{k_c}.'+eN{k_c}*yN{k_c}.';
        EUcorr{k_c} = EUcorr{k_c}-e1{k_c}*u1{k_c}.'+eN{k_c}*uN{k_c}.';
    end
    %}
    end

    % update plotting
    if rem(k2-1,update_interval)==0
        % top figure: update fill to indicate past data range
        if k2-1 >= update_interval
            delete(fillover);delete(xl1); delete(xl2);
            delete(im1u); delete(im2u); delete(im1y); delete(im2y);
        end
        xrange = [0.5+(k2-1) 0.5+Nbar+(k2-1)];
        fillover = fill(ax1,[xrange fliplr(xrange)],[y_max,y_max,y_min,y_min],fillcol,EdgeAlpha=0,FaceAlpha=0.15,HandleVisibility='off');
        ylim(ax1,y_lims);
        xl1 = xline(ax1,xrange(1),HandleVisibility='off');
        xl2 = xline(ax1,xrange(2),HandleVisibility='off');
        
        % EU correlation - regular DeePC
        im1u = imagesc(ax1u,EUcorr{1}/N{1}); colormap(ax1u,Cmap);colorbar(ax1u); xline(ax1u,p*nu+0.5,'k-',LineWidth=1.5);
        caxis(ax1u,[-cmaxu cmaxu]);
        ylim(ax1u,[0.5 0.5+ny*f]); set(ax1u, 'YDir','reverse');
        % EU correlation - CL-DeePC
        im2u = imagesc(ax2u,EUcorr{2}/N{2}); colormap(ax2u,Cmap);colorbar(ax2u); xline(ax2u,p*nu+0.5,'k-',LineWidth=1.5);
        caxis(ax2u,[-cmaxu cmaxu]);
        ylim(ax1y,[0.5 0.5+ny]);
        % EY correlation - regular DeePC
        im1y = imagesc(ax1y,EYcorr{1}/N{1}); colormap(ax1y,Cmap);colorbar(ax1y);
        caxis(ax1y,[-cmaxy cmaxy]);
        ylim(ax1y,[0.5 0.5+ny*f]); set(ax1y, 'YDir','reverse');
        % EY correlation - CL-DeePC
        im2y = imagesc(ax2y,EYcorr{2}/N{2}); colormap(ax2y,Cmap);colorbar(ax2y);
        caxis(ax2y,[-cmaxy cmaxy]);
        ylim(ax2y,[0.5 0.5+ny]);

        % pass figure to movie frame
        k3 = (k2-1)/update_interval+1;
        drawnow
        F(k3) = getframe(fig1);
%         F2 = frame2im(F2);
%         F2 = imresize(F2, 2);
%         F(k2) = im2frame(F2);
        writeVideo(vidfile,F(k3));
    end
end
close(vidfile)

%% helper functions
function [Efid,Yp,Upfid]=get_HankelMats(eCL,uCL,yCL,N,p,s,nu,ny) %s{k_c},N{k_c}
Nbar = length(eCL);

% initial hankel matrix - Ef
H_e  = mat2cell(eCL,ny,ones(1,Nbar));
H_e  = H_e(hankel(p+1:s,s:Nbar));
Efid = cell2mat(H_e);

% initial hankel matrix - Yp
H_y = mat2cell(yCL,ny,ones(1,Nbar));
H_y = H_y(hankel(1:p,p:p+N-1));
Yp  = cell2mat(H_y);

% initial hankel matrix - Upf
H_u = mat2cell(uCL,nu,ones(1,Nbar));
H_u = H_u(hankel(1:s,s:Nbar));
Upfid=cell2mat(H_u);
end

