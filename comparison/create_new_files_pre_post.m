%Rerun for each mouse to create new 'F_190128_GADB_plane_*_JA.mat' files for each plane 
%%
clear all

%Load selected_cells.mat
pre_post{1} = load('.\000-003\mGADd_000_000&mGADd_000_003_plane0_selected63.mat');
pre_post{2} = load('.\000-003\mGADd_000_000&mGADd_000_003_plane1_selected85.mat');
pre_post{3} = load('.\000-003\mGADd_000_000&mGADd_000_003_plane2_selected191.mat');
stim_post{1}= load('.\002-003\GADD_002&GADD_003_plane0_selected107.mat');
stim_post{2}= load('.\002-003\GADD_002&GADD_003_plane1_selected135.mat');
stim_post{3}= load('.\002-003\GADD_002&GADD_003_plane2_selected217.mat');
% change name for each mouse
post_matched = load('GADD_all_sessions_matched_cells.mat', 'matched_cells');
%%
%Load fall.mat
fall_000{1} = load('.\000\suite2p\plane0\Fall.mat');
fall_000{2} = load('.\000\suite2p\plane1\Fall.mat');
fall_000{3} = load('.\000\suite2p\plane2\Fall.mat');
fall_003{1} = load('.\003\suite2p\plane0\Fall.mat');
fall_003{2} = load('.\003\suite2p\plane1\Fall.mat');
fall_003{3} = load('.\003\suite2p\plane2\Fall.mat');

% change paths for each mouse
gad_mat_pre = load('.\000\mGADD_000_000.mat');
gad_mat_post = load('.\003\mGADD_000_003.mat');
ball_pre = load('.\000\mGADD_000_000_ball.mat', 'velocity');
ball_post = load('.\003\mGADD_000_003_ball.mat', 'velocity');
%%
%generate new file
d = dir('003\suite2p\plane*\Fall.mat');

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

    split_name = strsplit(d(p).folder, '\');
    recording_date = split_name{4};
    mouse_name = split_name{3};
    file_name = ['F_' recording_date '_' mouse_name '_plane_' num2str(p) '_JA.mat'];

%     %Pre and post stim 
    F_pre = fall_000{p}.F(pre_cell_num,:);
    Fneu_pre = fall_000{p}.Fneu(pre_cell_num,:);
    stat_pre = fall_000{p}.stat(pre_cell_num);
    ops_pre = fall_000{p}.ops;
    ops_pre.start_frame = gad_mat_pre.info.frame;
    ops_pre.mouse_name = mouse_name;
    velocity_pre = ball_pre.velocity;

    F_post = fall_003{p}.F(post_cell_num,:);
    Fneu_post = fall_003{p}.Fneu(post_cell_num,:);
    stat_post = fall_003{p}.stat(post_cell_num);
    ops_post = fall_003{p}.ops;
    ops_post.start_frame = gad_mat_post.info.frame;
    ops_post.mouse_name = mouse_name;
    velocity_post = ball_post.velocity;
    F_post_chan2 = fall_003{p}.F_post_chan2(post_cell_num, :);
    Fneu_post_chan2 = fall_003{p}.Fneu_post_chan2(post_cell_num, :);   

    for c = 1:size(F_pre, 1)
        stat_pre{c}.iplane = p;
        stat_post{c}.iplane = p;
    end

% both 
    redcell = fall_003{p}.redcell(post_cell_num, :);
    stim_cell_index = stim_cell_num;


%   for pre and post stim combined
  save(file_name, 'F_pre', 'Fneu_pre', 'stat_pre', 'ops_pre', 'F_post', 'Fneu_post', 'stat_post', 'ops_post', 'stim_cell_index', 'velocity_pre', 'velocity_post', 'redcell', 'F_post_chan2', 'Fneu_post_chan2', '-append');
    

end