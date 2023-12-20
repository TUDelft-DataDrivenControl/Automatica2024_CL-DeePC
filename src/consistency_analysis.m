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
num_c = 2;
% load an example to determine CL_sim_steps
fn = strcat('Nbar_',num2str(Nbar_all(1)),'.mat');
a = load(fn,'results');
CL_sim_steps = size(a.results.u_CL{1},2);

% iterate over noise variances Re
Quantiles = 0.1:0.2:0.9;
n_quant = length(Quantiles); idx_median = ceil(n_quant/2);%nan(n_quant,num_n)
eGu = cell(num_n,2,2,2);
% [eLu,eLy,eGu] = deal(cell(num_n,2,2,2)); % quantiles w/ median of cost
NormRef = nan(1,num_n);
for kn = 1:num_n
    Nbar = Nbar_all(kn);
    fn = strcat('Nbar_',num2str(Nbar),'.mat');
    
    % load results of all noise realizations for specific Re
    a = load(fn,'results');
    for k_dt = 2%1:2 % data type OL/adaptive CL
        switch k_dt
            case 1
                use_idx = 1;
            case 2
                use_idx = Nbar+1:CL_sim_steps;
        end
        for k_c = 1:2
            idx1 = sub2ind(size(eGu),kn,1,k_dt,k_c);
            idx2 = sub2ind(size(eGu),kn,2,k_dt,k_c);

%             eLu_run_mean = cellfun(@(x) mean(x(:,use_idx)),a.results.eLu(:,k_c));
%             eLu{idx1} = quantile(eLu_run_mean,Quantiles);
%             eLu{idx2} = mean(eLu_run_mean);
%             eLy_run_mean = cellfun(@(x) mean(x(:,use_idx)),a.results.eLy(:,k_c));
%             eLy{idx1} = quantile(eLy_run_mean,Quantiles);
%             eLy{idx2} = mean(eLy_run_mean);
            eGu_run_mean = cellfun(@(x) mean(x(:,use_idx)),a.results.eGu(:,k_c));
            eGu{idx1} = quantile(eGu_run_mean,Quantiles);
            eGu{idx2} = mean(eGu_run_mean);
        end
    end
end

%% Plotting

color = {[220,50,32]/256,[0,90,181]/256,[0.4660, 0.6740, 0.1880]};%{'#DC3220','#005AB5'};
figure('color','white','Visible','on');
% figure('Units', 'pixels', 'pos', [80 80 680 680],'color','white','Visible', 'on');
% TL = tiledlayout(fig6,1,1,'TileSpacing','tight');
% ax = cell(1,3);
FaceAlpha = 0.6/(idx_median-1);
for k_et = 3%1:3 % error type
    switch k_et
        case 1
            eLG = eLu;
            title_str = 'Error $\Gamma_f \mathcal{K}_p^\mathrm{u}$';
            ylab = '$||\Gamma_f \mathcal{K}_p^\mathrm{u}-\widehat{\Gamma_f \mathcal{K}_p^\mathrm{u}}||_\mathrm{F}^2$';
        case 2
            eLG = eLy;
            title_str = 'Error $\Gamma_f \mathcal{K}_p^\mathrm{y}$';
            ylab = '$||\Gamma_f \mathcal{K}_p^\mathrm{y}-\widehat{\Gamma_f \mathcal{K}_p^\mathrm{y}}||_\mathrm{F}^2$';
        case 3
            eLG = eGu;
            title_str = 'Error $\mathcal{T}_f^\mathrm{u}$';
            ylab = '$||\mathcal{T}_f^\mathrm{u}-\widehat{\mathcal{T}}_f^\mathrm{u}||_\mathrm{F}^2$';
    end
%     ax{k_et} = nexttile;
%     eLGq_OL1 = cell2mat(eLG(:,1,1,1)); % quantiles (OL -> OL data)
%     eLGq_OL2 = cell2mat(eLG(:,1,1,2));
    eLGq_CL1 = cell2mat(eLG(:,1,2,1));
    eLGq_CL2 = cell2mat(eLG(:,1,2,2));
%     eLGm_OL1 = cell2mat(eLG(:,2,1,1)); % means
%     eLGm_OL2 = cell2mat(eLG(:,2,1,2));
    eLGm_CL1 = cell2mat(eLG(:,2,2,1));
    eLGm_CL2 = cell2mat(eLG(:,2,2,2));
    for idx = 1:idx_median-1
        ff0 = fill([Nbar_all fliplr(Nbar_all)],[eLGq_CL1(:,idx).',  fliplr(eLGq_CL1(:,end-idx+1).')],  color{1},HandleVisibility='off');
        ff0.EdgeAlpha = 0.1; ff0.FaceAlpha = FaceAlpha; hold on;
        ff1 = fill([Nbar_all fliplr(Nbar_all)],[eLGq_CL2(:,idx).',  fliplr(eLGq_CL2(:,end-idx+1).')],  color{2},HandleVisibility='off');
        ff1.EdgeAlpha = 0.1; ff1.FaceAlpha = FaceAlpha;
    end
    % mean and median
    plot(Nbar_all,eLGq_CL1(:,idx_median), Color=color{1},LineWidth=1.5);
    plot(Nbar_all,eLGm_CL1,Color=color{1},  LineStyle='--',LineWidth=1.5);
    plot(Nbar_all,eLGq_CL2(:,idx_median), Color=color{2},LineWidth=1.5);
    plot(Nbar_all,eLGm_CL2,Color=color{2},  LineStyle='--',LineWidth=1.5);
    grid on;
    xlabel('$\bar{N}$','Interpreter','latex','FontSize',12);
    ylabel(ylab,'Interpreter','latex','FontSize',12);
    title(title_str,'Interpreter','latex','FontSize',14);
    legend({'DeePC (median)','DeePC (mean)','CL-DeePC (median)','CL-DeePC (mean)'},'Interpreter','latex','FontSize',12);
    h = gca; h.YScale = 'log';
end