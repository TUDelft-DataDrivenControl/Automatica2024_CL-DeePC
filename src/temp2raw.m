function temp2raw(num_e,num_c,run_dir,raw_dir,desc_runs,varargin)
results = struct;

%% save additional variables in results structure
nvarargin = length(varargin);
for k = 1:nvarargin
    var_name = inputname(nargin-nvarargin+k);
    results.(var_name) = varargin{k};
end

%% obtaining variables contained in directory description -> res2 structure
match = regexp(raw_dir,'Varying_([A-Za-z]+)_','tokens','once'); match = match{1};
dec_pat = '([-+]?\d*\.?\d+([eE][-+]?\d+)?)';
switch match
    case 'Nbar'
        str_pat = append('_p_(\d+)_f_(\d+)_Re_',dec_pat,'_Ru_',dec_pat,'_Rdu_',dec_pat,'_Q_',dec_pat,'_R_',dec_pat,'_dR_',dec_pat);
        matches = regexp(raw_dir,str_pat,'tokens','once');
        matches = cellfun(@str2double,matches,'UniformOutput',false);
        [res2.p,res2.f,res2.Re,res2.Ru,res2.Rdu,res2.Q,res2.R,res2.dR] = deal(matches{:});
    case 'pf'
        str_pat = append('_Nbar_(\d+)_Re_',dec_pat,'_Ru_',dec_pat,'_Rdu_',dec_pat,'_Q_',dec_pat,'_R_',dec_pat,'_dR_',dec_pat);
        matches = regexp(raw_dir,str_pat,'tokens','once');
        matches = cellfun(@str2double,matches,'UniformOutput',false);
        [res2.Nbar,res2.Re,res2.Ru,res2.Rdu,res2.Q,res2.R,res2.dR] = deal(matches{:});
    case 'Re'
        str_pat = append('_Nbar_(\d+)_p_(\d+)_f_(\d+)_Ru_',dec_pat,'_Rdu_',dec_pat,'_Q_',dec_pat,'_R_',dec_pat,'_dR_',dec_pat);
        matches = regexp(raw_dir,str_pat,'tokens','once');
        matches = cellfun(@str2double,matches,'UniformOutput',false);
        [res2.Nbar,res2.p,res2.f,res2.Ru,res2.Rdu,res2.Q,res2.R,res2.dR] = deal(matches{:});
    otherwise
        error('File naming not recognized');
end

%% saving variables in loaded results structure
% initialize results structure
results.CzLabel = {'DeePC, IV','CL-DeePC, IV','Oracle'}; % labels
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
results.eLu   = cell(num_e,2);
results.eLy   = cell(num_e,2);
results.eGu   = cell(num_e,2);
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
    
    fn1=sprintf('k%s_%d_ke_%d_ks_%d',letter,kvar,k_e,k_s);
    fn2=append(fn1,'.mat');
    temp_loc = fullfile(run_dir,fn2);
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
switch match
    case 'Nbar'
        load(temp_loc,'Nbar','ref');%,'ref','Ru','Rdu','Q','R','dR');
        results.Nbar = Nbar;
    case 'pf'
        load(temp_loc,'p','f','ref');%,'ref','Ru','Rdu','Q','R','dR');
        results.p = p;
        results.f = f;
    case 'Re'
        load(temp_loc,'Re','ref');%,'ref','Ru','Rdu','Q','R','dR');
        results.Re = Re;
end
results.ref = ref;

%[results.ref,results.Ru,results.Rdu,results.Q,results.R,results.dR,results.Rdu]=deal(ref,Ru,Rdu,Q,R,dR);

%% saving variables contained in directory description
fns2 = fieldnames(res2);
for k_fn2 = 1:length(fns2)
    fn2 = fns2{k_fn2};
    results.(fn2) = res2.(fn2);
end

%% save results structure
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