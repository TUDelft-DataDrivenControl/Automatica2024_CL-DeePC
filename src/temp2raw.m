function temp2raw(num_e,num_c,run_dir,raw_dir,desc_runs,varargin)

% save additional variables in results structure
nvarargin = length(varargin);
for k = 1:nvarargin
    var_name = inputname(nargin-nvarargin+k);
    results.(var_name) = varargin{k};
end

% initialize results structure
results.CzLabel = {'DeePC, IV','CL-DeePC, IV'}; % labels
results.ks    = nan(num_e,1);
results.noise = cell(num_e,1);       % innovation noise
results.du_CL = cell(num_e,1);
results.u_OL  = cell(num_e,1);
results.y_OL  = cell(num_e,1);
results.x_OL  = cell(num_e,1);
results.u_CL  = cell(num_e,num_c); % CL inputs
results.y_CL  = cell(num_e,num_c); % CL outputs
results.x_CL  = cell(num_e,num_c); % CL states
results.Cost  = cell(num_e,num_c);
results.eLu   = cell(num_e,num_c);
results.eLy   = cell(num_e,num_c);
results.eGu   = cell(num_e,num_c);
results.stat  = cell(num_e,num_c);

% Get a list of all .mat files in the directory
matFiles = dir(fullfile(run_dir, 'k*.mat'));

% Extract information using a single call to regexp
matches = regexp({matFiles.name}, 'k(.)_(\d+)_ke_(\d+)_ks_(\d+).mat', 'tokens', 'once');

% Convert cell array to matrix for numeric values
letter  = cellfun(@(x) x{1}, matches, 'UniformOutput', false); 
kvar = cellfun(@(x) str2double(x(2)), matches); 
ke_all = cellfun(@(x) str2double(x(3)), matches);
ks_all = cellfun(@(x) str2double(x(4)), matches);

% the below should be the same everywhere
letter   = letter{1};
kvar = kvar(1);

% Extract and sort the information
[~, sortingIndex] = sort(ke_all,2,"ascend");
ks_all = ks_all(sortingIndex);

% load and fill results in structure
for k_e = 1:num_e
    k_s = ks_all(k_e);
    results.ks(k_e,1) = k_s; % saving random seed number

    fn = sprintf('k%s_%d_ke_%d_ks_%d.mat',letter,kvar,k_e,k_s);
    temp_loc = fullfile(run_dir,fn);
    load(temp_loc,'e','du_CL','u_OL','y_OL','x_OL','u_CL','y_CL','x_CL','Cost','eLu','eLy','eGu','eObX','stat');
    results.noise{k_e,1} = e;
    results.du_CL{k_e,1} = du_CL;
    results.u_OL{k_e,1}  = u_OL;
    results.y_OL{k_e,1}  = y_OL;
    results.x_OL{k_e,1}  = x_OL;
    results.u_CL(k_e,:)  = u_CL; % CL inputs
    results.y_CL(k_e,:)  = y_CL; % CL outputs
    results.x_CL(k_e,:)  = x_CL; % CL states
    results.Cost(k_e,:)  = Cost;
    results.eLu(k_e,:)   = eLu;
    results.eLy(k_e,:)   = eLy;
    results.eGu(k_e,:)   = eGu;
    results.eObX(k_e,:)  = eObX;
    results.stat(k_e,:)  = stat;
end
% save results structure
save_filename = strcat(desc_runs,'.mat');
file_path = fullfile(raw_dir,save_filename);
save(file_path,"results")

% only delete temporary files when results have been saved
% for k_e = 1:num_e
%     k_s = ks_all(k_e);
%     fn = sprintf('k%s_%d_ke_%d_ks_%d.mat',letter,kvar,k_e,k_s);
%     temp_loc = fullfile(run_dir,fn);
%     delete(temp_loc);
% end
rmdir(run_dir,'s');

end