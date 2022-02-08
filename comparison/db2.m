function [plane1, plane2, plane3, epos] = db2()
% get cell data

plane1{1} = '.\GAD1\190214\F_190214_GAD1_plane_1_JA.mat';
plane2{1} = '.\GAD1\190214\F_190214_GAD1_plane_2_JA.mat';
plane3{1} = '.\GAD1\190214\F_190214_GAD1_plane_3_JA.mat';

plane1{2} = '.\GADA\190128\F_190128_GADA_plane_1_JA.mat';
plane2{2} = '.\GADA\190128\F_190128_GADA_plane_2_JA.mat';
plane3{2} = '.\GADA\190128\F_190128_GADA_plane_3_JA.mat';

plane1{3} = '.\GADB\190128\F_190128_GADB_plane_1_JA.mat';
plane2{3} = '.\GADB\190128\F_190128_GADB_plane_2_JA.mat';
plane3{3} = '.\GADB\190128\F_190128_GADB_plane_3_JA.mat';

plane1{4} = '.\GADC\190130\F_190128_GADD_plane_1_JA.mat';
plane2{4} = '.\GADC\190130\F_190128_GADD_plane_2_JA.mat';
plane3{4} = '.\GADC\190130\F_190128_GADD_plane_3_JA.mat';

plane1{5} = '.\GADD\190128\F_190130_GADC_plane_1_JA.mat';
plane2{5} = '.\GADD\190128\F_190130_GADC_plane_2_JA.mat';
plane3{5} = '.\GADD\190128\F_190130_GADC_plane_3_JA.mat';



% Load electrode positions
epos = load('.\epos.mat');
