%% get cell data
clear all;
d = dir('plane*/Fall.mat');
% ball_file = load('..\..\mGADA_000_002_ball.mat', 'velocity'); % change
% for each mouse
%% create new file 
for jj = 1:numel(d)  %for each plane
    
    clear stat
    d1 = load([d(jj).folder '/Fall.mat']);
    d.folder;
    file_name = strsplit(d(jj).folder, '\');
    mouse_name = file_name{3};
    exp_date = file_name{4};
    plane_name = ['_plane' num2str(jj)];
    filename= ['F_stim_' mouse_name '_' exp_date plane_name '_JA.mat'];

%     Fcell{1} = d1.F;
%     FcellNeu{1} = d1.Fneu;
%     ops = d1.ops;
%     velocity = decimate(ball_file.velocity, 3);
%     Fcell_chan2 = d1.Fcell_chan2;
%     FcellNeu_chan2 = d1.FcellNeu_chan2;
    
    for i = 1:size(d1.F, 1)
        if isfield(d1.stat{i},'neuropil_mask') 
            d1.stat{i} = rmfield(d1.stat{i}, 'neuropil_mask');
        end
        if isfield(d1.stat{i}, 'chan2_prob')
            d1.stat{i} = rmfield(d1.stat{i}, 'chan2_prob');
        end
        if isfield(d1.stat{i}, 'yext')
        else
            d1.stat{i}.yext = [];
        end
        if isfield(d1.stat{i}, 'xext')
        else    
            d1.stat{i}.xext = [];
        end

        stat(i) = d1.stat{i};
       
    end

    for m = 1:size(d1.F,1)
        stat(m).iscell = d1.iscell(m,1);
        stat(m).iplane = jj;
    end
    for n = 1:size(d1.redcell,1)
        stat(n).redcell = d1.redcell(n,1);
    end
    for o = size(d1.F,1)+1:size(d1.F,1)
        stat(n).redcell = 0;
    end
    
%     if isfield(ops, 'mimg1')
%         continue
%     else
%         ops.mimg1 = ops.meanImg;
%     end

    save(filename, 'stat', '-append')
%     save(filename,'Fcell','FcellNeu','ops','stat', 'velocity', 'Fcell_chan2', 'FcellNeu_chan2', '-append');
end
