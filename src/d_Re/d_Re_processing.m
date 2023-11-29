% Go to specific raw data set folder: data/raw/...
run_dir = pwd;

% Get a list of all contents in the directory
contents = dir(run_dir);

% Extract folder names
folderNames = {contents([contents.isdir]).name};

% Define the regular expression pattern
pattern = '^Re_(\d+(\.\d+)?)$';

% Apply regexp to the entire cell array
matches = regexp(folderNames, pattern, 'tokens', 'once');
matches = matches(~cellfun('isempty', matches)); % non-empty matches only

%%
num_n = numel(matches);
% results = repmat(struct('e',[],'u_OL',[],'du_CL',[],'y_OL',[],'x_OL',[],'u_CL',[],...
%     'y_CL',[],'x_CL',[],'Cost',[],'eLu',[],'eLy',[],'eGu',[],'eObX',[]),1,num_n);

% iterate over folders / noise variances Re
for kn = 1:num_n
    % data files have structure kn_%d_ke_%d_ks_%d.mat
    % -> find folder & files for specific kn (-> Re)
    [dir_Re,mat_fns,ke_all,ks_all,Re] = kn2dir(kn);

    % Initialize temporary structure array for this iteration
    tempData = struct();

    % iterate over files / noise realizations using parfor
    num_e = length(ke_all);
    tic
    parfor ce = 1:num_e
        ke = ke_all(ce);
        ks = ks_all(ce);
        mat_fn = mat_fns{ce};

        % load data
        s = load(fullfile(dir_Re, mat_fn));
        fns = fieldnames(s);
        for k_fld = 1:length(fns)
            fn = fns(k_fld); fn = fn{1};
            tempData(ce).(fn) = s.(fn);
        end
    end
    toc

    % Concatenate the temporary structure array to the main structure array
    results(kn).Re = Re;
    fns = fieldnames(tempData);
    for k_fld = 1:length(fns)
        fn = fns(k_fld); fn = fn{1};
        results(kn).(fn) = tempData.(fn);
    end
end
%%
function [folder,matFiles,ke,ks,Re] = kn2dir(kn)
    % Specify the path to the directory
    parentDirectory = pwd;
    
    % Define the regular expression pattern
    pattern = strcat('^kn_',num2str(kn),'_ke_(\d+)_ks_(\d+)\.mat');
    
    % Get a list of all .mat files in the directory
    matFiles = dir(fullfile(parentDirectory, '**/*.mat'));
    
    % Extract file names
    fileNames = {matFiles.name};
    
    % Apply regexp to the entire cell array
    matches = regexp(fileNames, pattern, 'tokens', 'once');
    
    % Find the indices of the non-empty matches
    indices = find(~cellfun('isempty', matches));
    
    % get folder
    folder = matFiles(indices(1)).folder;
    Re = regexp(folder,'Re_(\d+\.?\d*)$','tokens','once');
    Re = str2double(Re{1});

    % select all relevant .mat files
    matFiles = {matFiles(indices).name};


    ke = cellfun(@(x) str2double(x{1}),matches(indices));%,'UniformOutput',false);
    ks = cellfun(@(x) str2double(x{2}),matches(indices));%,'UniformOutput',false);
end

% function [Re,Cost] = get_Re_Cost(kn,ke)
%     
% end