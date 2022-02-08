clear all
d_files = dir('*plane*_JA.mat');

for i = 1:numel(d_files)
    stim_files(i) = load(fullfile(d_files(i).folder, d_files(i).name));
    for j = 1:size(stim_files(i).stat, 2)
        stim_files(i).stat(j).iplane = i
    end
end

for i = 1:numel(d_files)
    if isfield(stim_files(i).ops, 'mimg1')
        continue
    else
        stim_files(i).ops.mimg1 = stim_files(i).ops.meanImg
    end
end

% pre post files

% F_post = cat(1, stim_files(1).F_post, stim_files(2).F_post, stim_files(3).F_post);
% F_pre = cat(1, stim_files(1).F_pre, stim_files(2).F_pre, stim_files(3).F_pre);
% Fneu_post = cat(1, stim_files(1).Fneu_post, stim_files(2).Fneu_post, stim_files(3).Fneu_post);
% Fneu_pre = cat(1, stim_files(1).Fneu_pre, stim_files(2).Fneu_pre, stim_files(3).Fneu_pre);
% ops_pre{1} = stim_files(1).ops_pre;
% ops_pre{2} = stim_files(2).ops_pre;
% ops_pre{3} = stim_files(3).ops_pre;
% ops_post{1} = stim_files(1).ops_post;
% ops_post{2} = stim_files(2).ops_post;
% ops_post{3} = stim_files(3).ops_post;
% stat_post = [stim_files(1).stat_post, stim_files(2).stat_post, stim_files(3).stat_post];
% stat_pre = [stim_files(1).stat_pre, stim_files(2).stat_pre, stim_files(3).stat_pre];
% redcell = cat(1, stim_files(1).redcell, stim_files(2).redcell, stim_files(3).redcell);
% stim_cell_index = cat(2, stim_files(1).stim_cell_index, stim_files(2).stim_cell_index, stim_files(3).stim_cell_index);
% velocity_pre = stim_files(1).velocity_pre;
% velocity_post = stim_files(1).velocity_post;

% stim files

Fcell{1} = [stim_files(1).Fcell{1}; stim_files(2).Fcell{1}; stim_files(3).Fcell{1}];
FcellNeu{1} = [stim_files(1).FcellNeu{1}; stim_files(2).FcellNeu{1}; stim_files(3).FcellNeu{1}];
ops{1} = stim_files(1).ops;
ops{2} = stim_files(2).ops;
ops{3} = stim_files(3).ops;
stat = [stim_files(1).stat, stim_files(2).stat, stim_files(3).stat];

split_filename = strsplit(d_files(1).folder, '\');
mouse_name = split_filename{3};
file_name = ['F_stim_' mouse_name '_allpairs_JA.mat'];
save(file_name, 'Fcell', 'FcellNeu', 'ops', 'stat')
% save(file_name, 'F_post', 'F_pre', 'Fneu_post', 'Fneu_pre', 'ops_post', 'ops_pre', 'stat_post', 'stat_pre', 'redcell', 'stim_cell_index', 'velocity_post', 'velocity_pre');