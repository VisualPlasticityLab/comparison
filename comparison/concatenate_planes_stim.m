
clear all
d_files = dir('F_stim_*plane*_JA.mat');

for i = 1:numel(d_files)
    stim_files(i) = load(fullfile(d_files(i).folder, d_files(i).name),'stat');
end

% Fcell{1} = [stim_files(1).Fcell{1}; stim_files(2).Fcell{1}; stim_files(3).Fcell{1}];
% FcellNeu{1} = [stim_files(1).FcellNeu{1}; stim_files(2).FcellNeu{1}; stim_files(3).FcellNeu{1}];
% ops{1} = stim_files(1).ops;
% ops{2} = stim_files(2).ops;
% ops{3} = stim_files(3).ops;
stat = [stim_files(1).stat, stim_files(2).stat, stim_files(3).stat];
% Fcell_chan2 = [stim_files(1).Fcell_chan2;stim_files(2).Fcell_chan2;stim_files(3).Fcell_chan2];
% FcellNeu_chan2 = [stim_files(1).FcellNeu_chan2;stim_files(2).FcellNeu_chan2;stim_files(3).FcellNeu_chan2];
% velocity = stim_files(1).velocity;

split_filename = strsplit(d_files(1).folder, '\');
mouse_name = split_filename{3};
file_name = ['F_stim_' mouse_name '_allplanes_JA.mat'];

save(file_name, 'stat', '-append')
% save(file_name, 'Fcell', 'FcellNeu', 'ops', 'stat', 'Fcell_chan2', 'FcellNeu_chan2', 'velocity');