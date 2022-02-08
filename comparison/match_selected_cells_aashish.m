%Load relevant files

[pp0, folder1]=uigetfile('*.mat','select pre-post plane0 selected cells matfile');
[pp1, folder2]=uigetfile('*.mat','select pre-post plane1 selected cells matfile');
[pp2, folder3]=uigetfile('*.mat','select pre-post plane2 selected cells matfile');
[sp0, folder4]=uigetfile('*.mat','select stim-post plane0 selected cells matfile');
[sp1, folder5]=uigetfile('*.mat','select stim-post plane1 selected cells matfile');
[sp2, folder6]=uigetfile('*.mat','select stim-post plane2 selected cells matfile');

pre_post{1} = load(fullfile(folder1, pp0))
pre_post{2} = load(fullfile(folder2, pp1))
pre_post{3} = load(fullfile(folder3, pp2))
stim_post{1}= load(fullfile(folder4, sp0))
stim_post{2}= load(fullfile(folder5, sp1))
stim_post{3}= load(fullfile(folder6, sp2))

for m = 1:numel(pre_post)
    matched_cells{m} = [];
    tmp_matched_cells = [];
    for i = 1:numel(pre_post{m}.pair2)
        for j = 1:numel(stim_post{m}.pair2)
            if pre_post{m}.pair2(i) == stim_post{m}.pair2(j)
                tmp_matched_cells(end+1) = stim_post{m}.pair2(j);         
            end
        end
    matched_cells{m} = tmp_matched_cells;
    end
    fprintf('%d,%d,%d;',numel(pre_post{m}.pair2),numel(stim_post{m}.pair2),numel(matched_cells{m}))

end
% matched_cells_all{ii} = matched_cells

total_cells = numel(matched_cells{1})+ numel(matched_cells{2})+ numel(matched_cells{3}) 
section_name = strsplit(folder1, '\');
new_file_name = string([section_name(3), 'all_sessions_matched_cells.mat']);    
join_name = append(new_file_name{1},'_',new_file_name{2})
save(join_name, 'matched_cells', 'total_cells')