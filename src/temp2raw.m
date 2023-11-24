function temp2raw(num_e,num_c,temp_dir1,temp_dir2,raw_dir,desc_runs,varargin)

% save additional variables in results structure
nvarargin = length(varargin);
for k = 1:nvarargin
    var_name = inputname(nargin-nvarargin+k);
    results.(var_name) = varargin{k};
end

% initialize results structure
results.CzLabel = {'DeePC, IV','CL-DeePC, IV'}; % labels
results.noise = cell(1,num_e);       % innovation noise
results.du_CL = cell(1,num_e);
results.u_OL  = cell(1,num_e);
results.y_OL  = cell(1,num_e);
results.x_OL  = cell(1,num_e);
results.u_CL  = cell(1,num_e,num_c); % CL inputs
results.y_CL  = cell(1,num_e,num_c); % CL outputs
results.x_CL  = cell(1,num_e,num_c); % CL states
results.Cost  = cell(1,num_e,num_c);
results.eLu   = cell(1,num_e,num_c);
results.eLy   = cell(1,num_e,num_c);
results.eGu   = cell(1,num_e,num_c);

% load and fill results in structure
for k_e = 1:num_e
    temp_loc = fullfile(temp_dir2,strcat('ke_',num2str(k_e),'.mat'));
    load(temp_loc);
    results.noise{1,k_e} = e;
    results.du_CL{1,k_e} = du_CL;
    results.u_OL{1,k_e}  = u_OL;
    results.y_OL{1,k_e}  = y_OL;
    results.x_OL{1,k_e}  = x_OL;
    results.u_CL(1,k_e,:)  = u_CL; % CL inputs
    results.y_CL(1,k_e,:)  = y_CL; % CL outputs
    results.x_CL(1,k_e,:)  = x_CL; % CL states
    results.Cost(1,k_e,:)  = Cost;
    results.eLu(1,k_e,:)   = eLu;
    results.eLy(1,k_e,:)   = eLy;
    results.eGu(1,k_e,:)   = eGu;
    results.eObX(1,k_e,:)  = eObX;
    delete(temp_loc);
end
rmdir(temp_dir2);
rmdir(temp_dir1);

% save results structure
save_filename = strcat(desc_runs,'.mat');
file_path = fullfile(raw_dir,save_filename);
save(file_path,"results")

end