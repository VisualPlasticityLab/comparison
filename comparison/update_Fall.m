%%
%Author: Aashish Khimasia
%Date: 5th August 2021
clear all
% getting things from new (run later)
% new_fall{1} = '.\suite2p\plane0\Fall.mat';
% new_fall{2} = '.\suite2p\plane1\Fall.mat';
% new_fall{3} = '.\suite2p\plane2\Fall.mat';

% saving things to old 
old_fall{1} = '.\tiffs\suite2p\plane0\Fall.mat';
old_fall{2} = '.\tiffs\suite2p\plane1\Fall.mat';
old_fall{3} = '.\tiffs\suite2p\plane2\Fall.mat';

F_chan2_file{1} = readNPY('.\suite2p\plane0\F_chan2.npy');
F_chan2_file{2} = readNPY('.\suite2p\plane1\F_chan2.npy');
F_chan2_file{3} = readNPY('.\suite2p\plane2\F_chan2.npy');

Fneu_chan2_file{1} = readNPY('.\suite2p\plane0\Fneu_chan2.npy');
Fneu_chan2_file{2} = readNPY('.\suite2p\plane1\Fneu_chan2.npy');
Fneu_chan2_file{3} = readNPY('.\suite2p\plane2\Fneu_chan2.npy');

% redcell_file{1} = readNPY('.\suite2p\plane0\redcell.npy');
% redcell_file{2} = readNPY('.\suite2p\plane1\redcell.npy');
% redcell_file{3} = readNPY('.\suite2p\plane2\redcell.npy');
%% 
for i = 1:3
%     new = load(new_fall{i});
%     old = load(old_fall{i});

    % old.ops.meanImg_chan2 = new.ops.meanImg_chan2;
    % old.ops.meanImg_chan2_corrected = new.ops.meanImg_chan2;
%     redcell = redcell_file{i}; 

%     ops = old.ops;
%     stat = old.stat;
    Fcell_chan2 = F_chan2_file{i};
    FcellNeu_chan2 = Fneu_chan2_file{i};
    save(old_fall{i}, 'Fcell_chan2', 'FcellNeu_chan2', '-append')   
end