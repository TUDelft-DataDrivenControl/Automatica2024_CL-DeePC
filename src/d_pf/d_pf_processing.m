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
num_c = 3;
% iterate over noise variances Re
quantiles = 0.1:0.2:0.9;%[0.051 0.25 0.5 0.75 0.9];
n_quant = length(quantiles); idx_median = ceil(n_quant/2);
[costs_DeePC1,costs_CLDeePC1,costs_Oracle1] = deal(nan(n_quant,num_n)); % quantiles w/ median of cost
[costs_DeePC2,costs_CLDeePC2,costs_Oracle2] = deal(nan(1,num_n));       % mean cost

% pre-determine some constants
pf = pf_all(1);
fn = strcat('p_',num2str(pf),'_f_',num2str(pf),'.mat');
a = load(fn,'results');
len_x = length(a.results.y_CL{1});
NormRef = sqrt(mean(a.results.ref(:,Nbar+1:len_x).^2));
parfor kn = 1:num_n
    pf = pf_all(kn);
    fn = strcat('p_',num2str(pf),'_f_',num2str(pf),'.mat');
    
    % load results of all noise realizations for specific Re
    a = load(fn,'results');
%     costs = cellfun(@(x) mean(x(Nbar+1:end)),a.results.Cost);
    len_x = length(a.results.y_CL{1});
    costs = cellfun(@(x) sqrt(mean((x(Nbar+1:end)-a.results.ref(:,Nbar+1:len_x)).^2))/NormRef,a.results.y_CL);
    costs_DeePC1(:,kn)   = quantile(costs(:,1),quantiles);
    costs_CLDeePC1(:,kn) = quantile(costs(:,2),quantiles);
    costs_DeePC2(:,kn)   = mean(costs(:,1));
    costs_CLDeePC2(:,kn) = mean(costs(:,2));
    if num_c > 2
        costs_Oracle1(:,kn)  = quantile(costs(:,3),quantiles);
        costs_Oracle2(:,kn)  = mean(costs(:,3));
    end
end
%%
color = {[220,50,32]/256,[0,90,181]/256,[0.4660, 0.6740, 0.1880]};%{'#DC3220','#005AB5'};
fig3 = figure('Units', 'pixels', 'pos', [80 80 680 450],'color','white','Visible', 'on');
set(fig3,'Units','centimeters');
pos3 = get(fig3,'Position');
width3 = 8.4*1.5;
scale3 = width3/pos3(3);
set(fig3,'Position',[pos3(1:2),width3,scale3*pos3(4)])
[ff0,ff1,ff2] = deal(cell(1,idx_median-1));
FaceAlpha = 0.6/(idx_median-1);
for idx = 1:idx_median-1
    if num_c > 2
        ff0{idx}=fill([pf_all fliplr(pf_all)],[costs_Oracle1(idx,:),  fliplr(costs_Oracle1(end-idx+1,:))],color{3},HandleVisibility='off');
        ff0{idx}.EdgeAlpha = 0.1; ff0{idx}.FaceAlpha = FaceAlpha; hold on;
    end
    ff1{idx}=fill([pf_all fliplr(pf_all)],[costs_DeePC1(idx,:),  fliplr(costs_DeePC1(end-idx+1,:))],  color{1},HandleVisibility='off');
    ff1{idx}.EdgeAlpha = 0.1; ff1{idx}.FaceAlpha = FaceAlpha; hold on;
    ff2{idx}=fill([pf_all fliplr(pf_all)],[costs_CLDeePC1(idx,:),fliplr(costs_CLDeePC1(end-idx+1,:))],color{2},HandleVisibility='off');
    ff2{idx}.EdgeAlpha = 0.1; ff2{idx}.FaceAlpha = FaceAlpha; hold on;
end
if num_c > 2
    l0 = plot(pf_all,costs_Oracle1(idx_median,:),  Color=color{3},LineWidth=1.5); hold on;
    l01= plot(pf_all,costs_Oracle2,Color=color{3}, LineStyle='--',LineWidth=1.5);
end
l1 = plot(pf_all,costs_DeePC1(idx_median,:),   Color=color{1},LineWidth=1.5); hold on;
l11= plot(pf_all,costs_DeePC2,Color=color{1},  LineStyle='--',LineWidth=1.5);
l2 = plot(pf_all,costs_CLDeePC1(idx_median,:), Color=color{2},LineWidth=1.5);
l21= plot(pf_all,costs_CLDeePC2,Color=color{2},LineStyle='--',LineWidth=1.5);
% yline(sqrt(Re)./NormRef,'k-',LineWidth=1.5);
if num_c > 2
    leg = legend('Oracle (median)','Oracle (mean)','DeePC (median)','DeePC (mean)','CL-DeePC (median)','CL-DeePC (mean)');%,'ideal: $\sigma(e_k)/\sqrt{\overline{r_k^2}}$');
else
    leg = legend('DeePC (median)','DeePC (mean)','CL-DeePC (median)','CL-DeePC (mean)');%,'ideal: $\sigma(e_k)/\sqrt{\overline{r_k^2}}$');
end
set(leg,'Interpreter','latex','Location','northwest','ItemTokensize',13*ones(6,1),'FontSize',10); %use ItemTokensize to shorten legend lines
legpos = leg.Position;set(leg,'Position',[legpos(1:2) legpos(3)-0.02 legpos(4)]); % moves figure to the left a bit
grid on
% setting ticks based on y limits
% y_lim = ylim; 
% tickmarks = get_ticks(y_lim);
% yticks(tickmarks);
% yticklabels(tickmarks);
% ytickformat('%.0e');
% change_yticks(gca);
h =gca; h.YScale = 'log';

xlabel('$f(=p)$','Interpreter','latex','FontSize',13);
% ylabel('$\bar{J}$','Interpreter','latex','FontSize',12);
ylabel('$J_\mathrm{rms}$','Interpreter','latex','FontSize',13);
title(append('$\bar{N}=',num2str(Nbar),'$, ',...
      ...'$\sigma^2(d^\mathrm{u}_k)=',num2str(Rdu),'$, ',...
      '$\Sigma(e_k)=',num2str(Re),'$'),'Interpreter','latex','FontSize',13);
hAxes = gca;hAxes.TickLabelInterpreter = 'latex';

%% Helper functions

% get tick marks given y-axis limits
function tickmarks = get_ticks(y_lim)
y_min = y_lim(1); y_max = y_lim(2);
expmin = floor(log10(y_min));
y_min2 = round(ceil(y_min*10^(-expmin))/10^(-expmin),-expmin);
expmax = floor(log10(y_max));
y_max2 = round(floor(y_max*10^(-expmax))/10^(-expmax),-expmax);
nticks = (expmax-expmin)*10-(y_min2*10^(-expmin))+(y_max2*10^(-expmax));
tickmarks = nan(1,nticks);
tickmarks(1)   = y_min2;
tickmarks(end) = y_max2;
for k = 1:nticks-2
    step_size = 10^floor(log10(tickmarks(k)));
    % catches round-off imprecision error
    if tickmarks(k)/step_size >= 9.5
        step_size = step_size*10;
    end
    tickmarks(k+1) = tickmarks(k) + step_size;
end
end

% adjust y-axis tick labels
function change_yticks(ax)
data = ax.YTickLabels;

% get major tick labels
matches = regexp(data, '([\d.]+)\s*\\times10\^\{(-?\d+)\}', 'tokens', 'once');
mask1 = cellfun(@(x) (str2double(x{1})==5 || str2double(x{1})==1),matches); %str2double(x{1})==10 || 
mask2 = cellfun(@(x) str2double(x{1})==10, matches);
ax.YTick = ax.YTick((mask1|mask2));

data = ax.YTickLabels;
% Initialize a new cell array for modified data
modified_data = cell(size(data));

% Loop through each element in the original cell array
for i = 1:numel(data)
    % Extract the first number and exponent using regular expressions
    match = regexp(data{i}, '([\d.]+)\s*\\times10\^\{(-?\d+)\}', 'tokens', 'once');

    % Extracted values
    first_number = str2double(match{1});
    exponent = str2double(match{2});

    % Check if the first number is 10
    if first_number == 10
        exponent = exponent + 1;
        modified_data{i} = sprintf('10^{%d}', exponent);
    elseif first_number == 5
        modified_data{i} = sprintf('%d\\times10^{%d}', first_number, exponent);
    elseif first_number == 1
        modified_data{i} = sprintf('10^{%d}', exponent);
    else
        modified_data{i} = "";
    end
    
end
ax.YTickLabels = modified_data;
end