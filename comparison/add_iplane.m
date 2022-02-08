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

stat = [stim_files(1).stat, stim_files(2).stat, stim_files(3).stat];

ops{1} = stim_files(1).ops;
ops{2} = stim_files(2).ops;
ops{3} = stim_files(3).ops;

split_filename = strsplit(d_files(1).folder, '\');
mouse_name = split_filename{3};
file_name = ['F_stim_' mouse_name '_allpairs_JA.mat'];
save(file_name, 'stat', 'ops', '-append');