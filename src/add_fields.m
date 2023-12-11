clear
clc;

% go to directory with .mat files to add ref to
raw_dir = pwd;

% get all subdirectories that contain data
contents = dir(raw_dir);
contents = contents([contents.isdir]);
contents = contents(~strcmp({contents.name},'.' ));
contents = contents(~strcmp({contents.name},'..'));

for k_dir = 1:length(contents)
    sub_raw_dir = contents(k_dir).name;
    cd(sub_raw_dir);
    % Get a list of all .mat files in the directory
    matFiles = dir('*.mat');
    
    if ~isempty(matFiles)
        match = regexp(matFiles(1).folder,'Varying_([A-Za-z]+)_','tokens','once'); match = match{1};
        dec_pat = '([-+]?\d*\.?\d+([eE][-+]?\d+)?)';
        switch match
            case 'Nbar'
                str_pat = append('_p_(\d+)_f_(\d+)_Re_',dec_pat,'_Ru_',dec_pat,'_Rdu_',dec_pat,'_Q_',dec_pat,'_R_',dec_pat,'_dR_',dec_pat);
                matches = regexp(matFiles(1).folder,str_pat,'tokens','once');
                matches = cellfun(@str2double,matches,'UniformOutput',false);
                [res2.p,res2.f,res2.Re,res2.Ru,res2.Rdu,res2.Q,res2.R,res2.dR] = deal(matches{:});
            case 'pf'
                str_pat = append('_Nbar_(\d+)_Re_',dec_pat,'_Ru_',dec_pat,'_Rdu_',dec_pat,'_Q_',dec_pat,'_R_',dec_pat,'_dR_',dec_pat);
                matches = regexp(matFiles(1).folder,str_pat,'tokens','once');
                matches = cellfun(@str2double,matches,'UniformOutput',false);
                [res2.Nbar,res2.Re,res2.Ru,res2.Rdu,res2.Q,res2.R,res2.dR] = deal(matches{:});
            case 'Re'
                str_pat = append('_Nbar_(\d+)_p_(\d+)_f_(\d+)_Ru_',dec_pat,'_Rdu_',dec_pat,'_Q_',dec_pat,'_R_',dec_pat,'_dR_',dec_pat);
                matches = regexp(matFiles(1).folder,str_pat,'tokens','once');
                matches = cellfun(@str2double,matches,'UniformOutput',false);
                [res2.Nbar,res2.p,res2.f,res2.Ru,res2.Rdu,res2.Q,res2.R,res2.dR] = deal(matches{:});
            otherwise
                error('File naming not recognized');
        end

        parfor k_fn = 1:length(matFiles)
            save_new_results(matFiles,res2,k_fn)
        end
        clear res2
    end % of if

    cd('../');
end
function save_new_results(matFiles,res2,k_fn)
fn = matFiles(k_fn).name;
load(fn,'results');
fns2 = fieldnames(res2);
for k_fn2 = 1:length(fns2)
    fn2 = fns2{k_fn2};
    results.(fn2) = res2.(fn2);
end
save(fn,'results');
end