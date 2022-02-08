[file,path]=uigetfile('*signals.mat','select sigF file');
signals=load(fullfile(path,file))
matrix = signals.matrix;
win_sig = [signals.window(1) signals.window(end)];
sigF= signals.sigF;

[~,SI.peakR,SI.errorR,~,SI.peakS,SI.errorS]= sigFcmp(sigF,win_sig,matrix);  % sigR:seg,1,Var,ncell
%Cor= [14 16 19 4 115 37 28 27 120 34 116 25 118 117 23 32 43 33 9 10 13 49 47 61 59 112 123 57 55 54 41 86 85 87 88 90 77 75 74 72 76 79 64 63 96 99 101 107 92 108 110 68 69 70 71 108 29 8 11];
%Cor = [24 26 27 14 132 133 136 36 134 44 127 35 34 130 33 41 53 43 19 20 23 62 59 71 69 124 128 131 67 66 51 138 94 95 96 93 87 85 84 82 86 89 74 73 106 108 109 115 100 117 120 79 80 81 129 117 37 18 21]
%Cor = [ 2 13 20 21 29 62 6 9 17 18 23 60 65 27 28 61 63 31 30 40 53 67 55 42 34 32 57 58 59 51 52 44 45 36 37 35];
Cor = [ 72 71 22 23 38 30 5 8 15 16 25 26 27 35 70 29 19 69 39 46 61 73 63 49 41 40 65 66 67 59 60 52 53 43 44 42];
fixpeak(SI.peakR(:,:,Cor));