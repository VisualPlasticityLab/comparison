% Author = Aashish Khimasia
% Date = 06/08/21
%% Load epos and matched cell data across all planes for each mouse
function [epos, se] = create_se()
epos = load('epos.mat');
se{1} = load('F_GAD1_allplanes_JA.mat');
se{2} = load('F_GADA_allplanes_JA.mat');
se{3} = load('F_GADB_allplanes_JA.mat');
se{4} = load('F_GADC_allplanes_JA.mat');
se{5} = load('F_GADD_allplanes_JA.mat');