% go to directory with folders to convert to data
raw_dir = pwd;

% get all subdirectories that contain data
contents = dir(raw_dir);
contents = contents([contents.isdir]);
contents = contents(~strcmp({contents.name},'.' ));
contents = contents(~strcmp({contents.name},'..'));

% determine num_e & num_c by looking at one data file
if ~isempty(contents)
    matFiles = dir(fullfile(raw_dir,contents(1).name,'*.mat'));
    load(fullfile(matFiles(1).folder,matFiles(1).name),'u_CL');
    num_e = numel(matFiles);
    num_c = numel(u_CL);
end

% run tem2raw on contents of each folder
for k_dir = 1:length(contents)
    run_dir = contents(k_dir).name;
    desc_runs = run_dir; % save .mat file with same name as folder
    
    % save all data in folder in a results structure and delete folder
    temp2raw(num_e,num_c,run_dir,raw_dir,desc_runs);
end