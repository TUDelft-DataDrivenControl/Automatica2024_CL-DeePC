%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Closed-loop Data-Enabled Predictive Control
%           Authors: R. Dinkla, S.P. Mulders, T. Oomen, J.W. van Wingerden
%           Submission to Automatica 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NB: first move into the main project directory
close all;
clear;
rng('default');
clc;

dec_pat = '([-+]?\d*\.?\d+([eE][-+]?\d+)?)';
main_dir = pwd;

%% Figure 1: J vs Re - adjust files used
clearvars -except main_dir dec_pat
dirname = 'Varying_Re_0.0001-1-50_Nbar_239_p_20_f_20_Ru_1_Rdu_0_Q_100_R_0_dR_10';
cd(fullfile('data','raw',dirname));
d_Re_processing;
fig1 = gcf;
set(fig1,'Color','w');
cd(fullfile('..','..','..','results','figures'));
exportgraphics(fig1,append(dirname,'.pdf'),'ContentType','vector');
cd(main_dir);

%% Figure 2: J vs Nbar
clearvars -except main_dir dec_pat
dirname = 'Varying_Nbar_99-1039-50_p_20_f_20_Re_1_Ru_1_Rdu_0_Q_100_R_0_dR_10';
% dirname = 'Varying_Nbar_99-1039-50_p_20_f_20_Re_0.25_Ru_1_Rdu_0_Q_100_R_0_dR_10';
cd(fullfile('data','raw',dirname));
d_Nbar_processing;
fig2 = gcf;
set(fig2,'Color','w');
cd(fullfile('..','..','..','results','figures'))
exportgraphics(fig2,append(dirname,'.pdf'),'ContentType','vector');
cd(main_dir);

%% Figure 3: J vs p=f
clearvars -except main_dir dec_pat
dirname = 'Varying_pf_20-100-41_Nbar_1199_Re_0.25_Ru_1_Rdu_0_Q_100_R_0_dR_10';
cd(fullfile('data','raw',dirname));
d_pf_processing;
fig3 = gcf;
set(fig3,'Color','w');
cd(fullfile('..','..','..','results','figures'))
exportgraphics(fig3,append(dirname,'.pdf'),'ContentType','vector');
cd(main_dir);

%% Fig 4: problem with regular DeePC
clearvars -except main_dir dec_pat
dirname = 'Varying_Re_0.0001-1-50_Nbar_239_p_20_f_20_Ru_1_Rdu_0_Q_100_R_0_dR_10';
fn      = 'Re_1.mat';
load(fullfile('data','raw',dirname,fn),'results');
max_costs = cellfun(@max,results.Cost(:,1));
[~,idx_sorted] = sort(max_costs,'descend');
k_e = idx_sorted(1);
color = {[220,50,32]/256,[0,90,181]/256,[0.4660, 0.6740, 0.1880]};

fig4 = figure('Units', 'pixels', 'pos', [80 80 680 450],'color','white','Visible', 'on');
ax1 = axes(fig4);
set(fig4,'Units','centimeters');
pos4 = get(fig4,'Position');
width4 = 8.4*1.5;
scale4 = width4/pos4(3);
set(fig4,'Position',[pos4(1:2),width4,scale4*pos4(4)])
plot(ax1,results.y_CL{k_e,1}.','color',color{1}); hold on;
plot(ax1,results.y_CL{k_e,2}.','color',color{2});
plot(ax1,results.y_CL{k_e,3}.','color',color{3});
plot(ax1,results.ref,'k');
legend(ax1,'DeePC','CL-DeePC','oracle','reference','Interpreter','latex','FontSize',10);
ylabel(ax1,'$y_k$','Interpreter','latex','FontSize',13);
grid on;
xlim(ax1,[0 length(results.ref)]);
ylim(ax1,[-70 450]);
xline(ax1,results.Nbar+0.5,HandleVisibility='off')
xlabel(ax1,'Time index $k$','Interpreter','latex','FontSize',13);
title_str = append('$\bar{N}=',num2str(results.Nbar),'$, ',...
      '$\Sigma(e_k)=',num2str(results.Re),'$, ',...
      ...'$\Sigma(d^\mathrm{u}_k)=',num2str(results.Rdu),'$, ',...
      '$p=',num2str(results.p),'$, $f=',num2str(results.f),'$');
title(title_str,'Interpreter','latex','FontSize',13);
hAxes = gca;hAxes.TickLabelInterpreter = 'latex';
set(fig4,'Color','w');

desc_str1 = regexp(fn,'(.*)\.mat$','tokens','once'); desc_str1 = desc_str1{1};
desc_str2 = regexp(dirname,'Nbar_(.*)$','tokens','once'); desc_str2 = desc_str2{1};
figname = append('DeePC_CL_ID_issue_',desc_str1,'_Nbar_',desc_str2,'.pdf');
cd(fullfile('results','figures'));
exportgraphics(fig4,figname,'ContentType','vector');
cd(main_dir);
%% Fig 5: adaptive CL-correlation between noise & inputs
clearvars -except main_dir dec_pat
dirname = 'Varying_Nbar_99-1039-50_p_20_f_20_Re_1_Ru_1_Rdu_0_Q_100_R_0_dR_10';
fn      = 'Nbar_1039.mat';
load(fullfile('data','raw',dirname,fn),'results');
corr_analysis2;
fig5=gcf;
set(fig5,'Color','w');
cd(fullfile('results','figures'));
matches = regexp(dirname,'Varying_Nbar_(\d+-\d+-\d+)_p_(.*)','tokens','once');
figname = append('Correlation_Nbar_',num2str(Nbar),'_p_',matches{2},'.pdf');
exportgraphics(fig5,figname,'ContentType','vector');
cd(main_dir);

%% Fig 6: Consistency analysis in adaptive CL applications
clearvars -except main_dir dec_pat
dirname = 'Varying_Nbar_99-1039-50_p_20_f_20_Re_1_Ru_1_Rdu_0_Q_100_R_0_dR_10';
cd(fullfile('data','raw',dirname));
consistency_analysis;
fig6=gcf;
set(fig6,'Color','w');
cd(fullfile('..','..','..','results','figures'))
matches = regexp(dirname,'Varying_Nbar_(\d+-\d+-\d+)_p_(.*)','tokens','once');
figname = append('Consistency_Nbar_',matches{1},'_p_',matches{2},'.pdf');
exportgraphics(fig6,figname,'ContentType','vector');
cd(main_dir);