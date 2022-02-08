clear all;close all;clc
% func_save_data2('day2')
% func_save_data2('day1')

figname=uigetfile('*.fig','find fig with the labeled cells');
open(figname);

figname2=uigetfile('*.fig','find fig with the labeled cells');
open(figname2);
%%
[day1eye,path]=uigetfile('*signals.mat','select trace for day1eye');
day1=load(fullfile(path,day1eye))
[day1eye2,path]=uigetfile('peakSI.mat','select peakSI matfile for day1eye');
resp1=load(fullfile(path,day2eye2));




[day2eye,path]=uigetfile('*signals.mat','open trace for day2eye');
day2=load(fullfile(path,day2eye))
[day2eye2,path]=uigetfile('peakSI.mat','select peakSI matfile for day2eye');
resp2=load(fullfile(path,day2eye2));


%% 
ncells = size(day1.sigF,1);
cellnum=1:ncells-1;
selected=inputdlg('input good cell#','includes:',[4 80]);
goodcell=str2num(selected{1});
goodcell1= goodcell(1:2:end);
goodcell2= goodcell(2:2:end);
[m,n] = size(goodcell1);
nr_gd = n;

ans1 = resp1.SI.peakR(:,:,goodcell1);
ans2 = resp2.SI.peakR(:,:,goodcell2);

for i = 1:nr_gd
    if any(ans1(:,:,i)) || any(ans2(:,:,i))
        [h(i),p(i)] = ttest(ans1(:,:,i),ans2(:,:,i),'Alpha',0.05)
        final_nr1(i) = goodcell1(i)
        final1(i,:) = ans1(:,:,i)
        final_nr2(i) = goodcell2(i)
        final2(i,:) = ans2(:,:,i)
    end
end
