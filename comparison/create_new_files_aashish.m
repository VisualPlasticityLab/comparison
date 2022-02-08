%Rerun for each mouse to create new mat files for each plane
%%
clear all
%Get selected_cells.mat
[pp0, folder1]=uigetfile('*.mat','select pre-post plane0 selected cells matfile');
[pp1, folder2]=uigetfile('*.mat','select pre-post plane1 selected cells matfile');
[pp2, folder3]=uigetfile('*.mat','select pre-post plane2 selected cells matfile');
[sp0, folder4]=uigetfile('*.mat','select stim-post plane0 selected cells matfile');
[sp1, folder5]=uigetfile('*.mat','select stim-post plane1 selected cells matfile');
[sp2, folder6]=uigetfile('*.mat','select stim-post plane2 selected cells matfile');

%%
%Load selected_cells.mat
pre_post{1} = load(fullfile(folder1, pp0));
pre_post{2} = load(fullfile(folder2, pp1));
pre_post{3} = load(fullfile(folder3, pp2));
stim_post{1}= load(fullfile(folder4, sp0));
stim_post{2}= load(fullfile(folder5, sp1));
stim_post{3}= load(fullfile(folder6, sp2));

post_matched = load('GAD1_all_sessions_matched_cells.mat', 'matched_cells');
% [GADmat, folder7]=uigetfile('*.mat','select GAD*_all_sessions_matched_cells matfile');
% post_matched = load(fullfile(folder7, GADmat));
%%
%Load fall.mat
fall_000{1} = load('Z:\Teams-MJ-data\GADB\190128\002\suite2p\plane0\Fall.mat');
fall_000{2} = load('Z:\Teams-MJ-data\GADB\190128\002\suite2p\plane1\Fall.mat');
fall_000{3} = load('Z:\Teams-MJ-data\GADB\190128\002\suite2p\plane2\Fall.mat');
fall_003{1} = load('Z:\Teams-MJ-data\GADB\190128\005\suite2p\plane0\Fall.mat');
fall_003{2} = load('Z:\Teams-MJ-data\GADB\190128\005\suite2p\plane1\Fall.mat');
fall_003{3} = load('Z:\Teams-MJ-data\GADB\190128\005\suite2p\plane2\Fall.mat');

gad_mat_pre = load('Z:\Teams-MJ-data\GADB\190128\002\mGADb_000_002.mat');
gad_mat_post = load('Z:\Teams-MJ-data\GADB\190128\005\mGADb_000_005.mat');

%%
%generate new file
post_cell_num = [];
stim_cell_num = [];
pre_cell_num = [];

for p = 1:numel(fall_000)
    post_cell_num = [];
    stim_cell_num = [];
    pre_cell_num = [];
    
    post_cells = post_matched.matched_cells{p};
    for j = 1:numel(post_cells)
        post_cell_num(end+1) = post_cells(j);
        stim_cell_num(end+1) = stim_post{p}.pair1(stim_post{p}.pair2==post_cells(j));
        pre_cell_num(end+1) = pre_post{p}.pair1(pre_post{p}.pair2==post_cells(j));
    end

    split_name = strsplit(folder1, '\');
    recording_date = split_name{4};
    mouse_name = split_name{3};
    file_name = ['F_' recording_date '_' mouse_name '_plane_' num2str(p) '_JA.mat'];

%     %Pre and post stim 
%     F_pre = fall_000{p}.F(pre_cell_num,:);
%     Fneu_pre = fall_000{p}.Fneu(pre_cell_num,:);
%     stat_pre = fall_000{p}.stat(pre_cell_num);
%     ops_pre = fall_000{p}.ops;
%     ops_pre.start_frame = gad_mat_pre.info.frame;
%     ops_pre.mouse_name = mouse_name;
%     ops_pre.plane_num = num2str(p);
%     
%     F_post = fall_003{p}.F(post_cell_num,:);
%     Fneu_post = fall_003{p}.Fneu(post_cell_num,:);
%     stat_post = fall_003{p}.stat(post_cell_num);
%     ops_post = fall_003{p}.ops;
%     ops_post.start_frame = gad_mat_post.info.frame;
%     ops_post.mouse_name = mouse_name;
%     ops_post.plane_num = num2str(p);

%   for pre and post stim combined
%   save(file_name, 'F_pre', 'Fneu_pre', 'stat_pre', 'ops_pre', 'F_post', 'Fneu_post', 'stat_post', 'ops_post');
    
%   for stim 

end
%% previous version
% post_cell_num = []
% stim_cell_num = []
% pre_cell_num = []
% %plane num + 1
% p = 1
% post_cells = post_matched.matched_cells{p}
% for i = 1:numel(post_cells)
% post_cell_num(end+1) = post_cells(i);
% stim_cell_num(end+1) = stim_post{p}.pair1(stim_post{p}.pair2==post_cells(i));
% pre_cell_num(end+1) = pre_post{p}.pair1(pre_post{p}.pair2==post_cells(i));
% end

