% Go to specific raw data set folder: data/raw/...
run_dir = pwd;

% Get a list of all contents in the directory
contents = dir(run_dir);

% Extract folder names
fileNames = {contents(~[contents.isdir]).name};

% Define the regular expression pattern
pattern = '^p_(\d+)_f_\d+.mat$';

% Apply regexp to the entire cell array
matches = regexp(fileNames, pattern, 'tokens', 'once');
matches = matches(~cellfun('isempty', matches)); % non-empty matches only
pf_all = cellfun(@(x) str2double(x{1}),matches);
pf_all = sort(pf_all,2,'ascend');

% get Nbar
Nbar = regexp(run_dir,'.*_Nbar_(\d+)_.*','tokens','once');
Nbar = str2double(Nbar{1});
% get Re
Re = regexp(run_dir,'.*_Re_(\d+(\.\d+)?)_.*','tokens','once');
Re = str2double(Re{1});
% get Ru
Ru = regexp(run_dir,'.*_Ru_(\d+(\.\d+)?)_.*','tokens','once');
Ru = str2double(Ru{1});
% get dRu
Rdu = regexp(run_dir,'.*_Rdu_(\d+(\.\d+)?)_.*','tokens','once');
Rdu = str2double(Rdu{1});
% get Q
Q = regexp(run_dir,'.*_Q_(\d+(\.\d+)?)_.*','tokens','once');
Q = str2double(Q{1});
% get dR
dR = regexp(run_dir,'.*_dR_(\d+(\.\d+)?).*','tokens','once');
dR = str2double(dR{1});

%%
num_n = numel(matches);

% iterate over noise variances Re
quantiles = 0.1:0.1:0.9;%[0.1 0.25 0.5 0.75 0.9];
n_quant = length(quantiles); idx_median = ceil(n_quant/2);
costs_DeePC1 = nan(n_quant,num_n); % quantiles w/ median of cost
costs_CLDeePC1 = costs_DeePC1;
costs_DeePC2 = nan(1,num_n);       % mean cost
costs_CLDeePC2 = costs_DeePC2;
parfor kn = 1:num_n
    pf = pf_all(kn);
    fn = strcat('p_',num2str(pf),'_f_',num2str(pf),'.mat');
    
    % load results of all noise realizations for specific Re
    a = load(fn,'results');
    costs = cellfun(@(x) mean(x(Nbar+1:end)),a.results.Cost);
    costs_DeePC1(:,kn)   = quantile(costs(:,1),quantiles);
    costs_CLDeePC1(:,kn) = quantile(costs(:,2),quantiles);
    costs_DeePC2(:,kn)   = mean(costs(:,1));
    costs_CLDeePC2(:,kn) = mean(costs(:,2));
end
%%
color = {[220,50,32]/256,[0,90,181]/256};%{'#DC3220','#005AB5'};
figure();
l1 = semilogy(pf_all,costs_DeePC1(idx_median,:),   Color=color{1},LineWidth=2.5); hold on;
l11= semilogy(pf_all,costs_DeePC2,Color=color{1},  LineStyle='--',LineWidth=1.5);
l2 = semilogy(pf_all,costs_CLDeePC1(idx_median,:), Color=color{2},LineWidth=2.5);
l21= semilogy(pf_all,costs_CLDeePC2,Color=color{2},LineStyle='--',LineWidth=1.5);
leg = legend('DeePC (median)','DeePC (mean)','CL-DeePC (median)','CL-DeePC (mean)');
set(leg,'Interpreter','latex');
ff1 = cell(1,idx_median-1);
ff2 = ff1;
FaceAlpha = 0.6/(idx_median-1);
for idx = 1:idx_median-1
    ff1{idx}=fill([pf_all fliplr(pf_all)],[costs_DeePC1(idx,:),  fliplr(costs_DeePC1(end-idx+1,:))],  Color{1},HandleVisibility='off');
    ff1{idx}.EdgeAlpha = 0.1; ff1{idx}.FaceAlpha = FaceAlpha; hold on;
    ff2{idx}=fill([pf_all fliplr(pf_all)],[costs_CLDeePC1(idx,:),fliplr(costs_CLDeePC1(end-idx+1,:))],Color{2},HandleVisibility='off');
    ff2{idx}.EdgeAlpha = 0.1; ff2{idx}.FaceAlpha = FaceAlpha; hold on;
end
grid on
xlabel('$f(=p)$','Interpreter','latex');
ylabel('$\bar{J}$','Interpreter','latex');
title(append('$\bar{N}=',num2str(Nbar),'$, ',...
      '$\sigma^2(e_k)=',num2str(Re),'$, ',...
      '$\sigma^2(d^\mathrm{u}_k)=',num2str(Rdu),'$'),'Interpreter','latex');
