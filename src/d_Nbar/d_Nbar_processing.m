% Go to specific raw data set folder: data/raw/...
run_dir = pwd;

% Get a list of all contents in the directory
contents = dir(run_dir);

% Extract folder names
fileNames = {contents(~[contents.isdir]).name};

% Define the regular expression pattern
pattern = '^Nbar_(.*).mat$';

% Apply regexp to the entire cell array
matches = regexp(fileNames, pattern, 'tokens', 'once');
matches = matches(~cellfun('isempty', matches)); % non-empty matches only
Nbar_all = cellfun(@(x) str2double(x{1}),matches);
Nbar_all = sort(Nbar_all,2,'ascend');
Nbar_max = Nbar_all(end);

% get Re
Re = regexp(run_dir,'.*_Re_(\d+(\.\d+)?)_.*','tokens','once');
Re = str2double(Re{1});
% get p
p = regexp(run_dir,'.*_p_(\d+)_.*','tokens','once');
p = str2double(p{1});
% get f
f = regexp(run_dir,'.*_f_(\d+)_.*','tokens','once');
f = str2double(f{1});
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
n_quant = length(quantiles); idx_mean = ceil(n_quant/2);
costs_DeePC1 = nan(n_quant,num_n); % quantiles w/ median of cost
costs_CLDeePC1 = costs_DeePC1;
costs_DeePC2 = nan(1,num_n);       % mean cost
costs_CLDeePC2 = costs_DeePC2;
parfor kn = 1:num_n
    Nbar = Nbar_all(kn);
    fn = strcat('Nbar_',num2str(Nbar),'.mat');
    
    % load results of all noise realizations for specific Re
    a = load(fn,'results');
    costs = cellfun(@(x) mean(x(Nbar+1:end)),a.results.Cost); % Nbar changes -> diff. #samples for mean! account for larger variance of mean?
    costs_DeePC1(:,kn)   = quantile(costs(:,1),quantiles);
    costs_CLDeePC1(:,kn) = quantile(costs(:,2),quantiles);
    costs_DeePC2(:,kn)   = mean(costs(:,1));
    costs_CLDeePC2(:,kn) = mean(costs(:,2));
end
%%
color = {[220,50,32]/256,[0,90,181]/256};%{'#DC3220','#005AB5'};
figure();
l1 = semilogx(Nbar_all,costs_DeePC1(idx_median,:),   Color=color{1},LineWidth=2.5); hold on;
l11= semilogx(Nbar_all,costs_DeePC2,Color=color{1},  LineStyle='--',LineWidth=1.5);
l2 = semilogx(Nbar_all,costs_CLDeePC1(idx_median,:), Color=color{2},LineWidth=2.5);
l21= semilogx(Nbar_all,costs_CLDeePC2,Color=color{2},LineStyle='--',LineWidth=1.5);
leg = legend('DeePC (median)','DeePC (mean)','CL-DeePC (median)','CL-DeePC (mean)');
set(leg,'Interpreter','latex');
ff1 = cell(1,idx_mean-1);
ff2 = ff1;
FaceAlpha = 0.6/(idx_mean-1);
for idx = 1:idx_median-1
    ff1{idx}=fill([Nbar_all fliplr(Nbar_all)],[costs_DeePC1(idx,:),  fliplr(costs_DeePC1(end-idx+1,:))],  color{1},HandleVisibility='off');
    ff1{idx}.EdgeAlpha = 0.1; ff1{idx}.FaceAlpha = FaceAlpha; hold on;
    ff2{idx}=fill([Nbar_all fliplr(Nbar_all)],[costs_CLDeePC1(idx,:),fliplr(costs_CLDeePC1(end-idx+1,:))],color{2},HandleVisibility='off');
    ff2{idx}.EdgeAlpha = 0.1; ff2{idx}.FaceAlpha = FaceAlpha; hold on;
end
grid on
xlim([Nbar_all(1) Nbar_all(end)])
xlabel('$\bar{N}$','Interpreter','latex');
ylabel('$\bar{J}$','Interpreter','latex');
title(append('$\sigma^2(e_k)=',num2str(Re),'$, ',...
      '$\sigma^2(d^\mathrm{u}_k)=',num2str(Rdu),'$, ',...
      '$p=',num2str(p),'$, $f=',num2str(f),'$'),'Interpreter','latex');
