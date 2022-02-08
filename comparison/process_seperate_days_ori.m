%% Load Day1 and Day2 traces and peakSI info
% [day1sigf,path]=uigetfile('*signals.mat','select trace for day1eye');
% day1sig=load(fullfile(path,day1sigf))
% day1.sigF= permute(day1sig.sigF,[4 3 2 1]);
[day1peakf,path1]=uigetfile('peakSI.mat','select peakSI matfile for day1eye');
[day2peakf,path2]=uigetfile('peakSI.mat','select peakSI matfile for day2eye');

%%
day1peak=load(fullfile(path1,day1peakf));
day1.finalvalue = permute(day1peak.SI.peakS,[3 2 1]);
day1.stdofeachvalue = permute(day1peak.SI.errorS,[3 2 1]);
day1.finalvalue2 = permute(day1peak.SI.peakR,[3 2 1]);
day1.stdofeachvalue2 = permute(day1peak.SI.errorR,[3 2 1]);
try
    day1.sigF= permute(day1peak.sigF,[4 3 2 1]);
catch
    day1.sigF= permute(day1peak.sigspF,[4 3 2 1]);
end
% [day2sigf,path]=uigetfile('*signals.mat','open trace for day2eye');
% day2sig=load(fullfile(path,day2sigf));
% day2.sigF= permute(day2sig.sigF,[4 3 2 1]);
day2peak=load(fullfile(path2,day2peakf));
day2.finalvalue = permute(day2peak.SI.peakS,[3 2 1]);
day2.stdofeachvalue = permute(day2peak.SI.errorS,[3 2 1]);
day2.finalvalue2 = permute(day2peak.SI.peakR,[3 2 1]);
day2.stdofeachvalue2 = permute(day2peak.SI.errorR,[3 2 1]);
try
    day2.sigF= permute(day2peak.sigF,[4 3 2 1]);
    suffix = '_CA';
catch
    day2.sigF= permute(day2peak.sigspF,[4 3 2 1]);
    suffix = '_SP';
end
%% sigF:ncells,nSteps, numberofcycles,framestocapture

nSteps=size(day1.sigF,2);
framestocapture=min (size(day1.sigF,4), size(day2.sigF,4));
nSteps1=nSteps-mod(nSteps,2);

[filename,path]=uigetfile('*.mat','Pick the preselected good cell?');
if filename
    pre = load(fullfile(path,filename));
    try
        goodcell1 = pre.selectedcell(:,1);
        goodcell2 = pre.selectedcell(:,2);
    catch
        goodcell1 = pre.pair1;
        goodcell2 = pre.pair2;
    end
else
    goodcell1 = 1:size(day1.finalvalue,1);
    goodcell2 = 1:size(day2.finalvalue,1);
end
%%
peak_selected1 =[];
peak_selected2 =[];

for j=1:numel(goodcell1)
    fit.peak1S = day1.finalvalue(goodcell1(j),:);
    fit.peak1R = day1.finalvalue2(goodcell1(j),:);
    peak_selected1(:,:,j)=[fit.peak1S;fit.peak1R];
end
for j=1:numel(goodcell2)
    fit.peak2S = day2.finalvalue(goodcell2(j),:);
    fit.peak2R = day2.finalvalue2(goodcell2(j),:);
    peak_selected2(:,:,j)=[fit.peak2S;fit.peak2R];
end
%%
if numel(goodcell1)~= numel(goodcell2)
    ncell = min(numel(goodcell1),numel(goodcell2));
    picked = unique(randi(ncell,1,ncell));
    selectedcell = [picked picked];
    peak_selected = [peak_selected1(:,:,picked);peak_selected2(:,:,picked)];
    h0 = plotpeakselected_new(peak_selected)
    figname=sprintf('random_%dcells',numel(picked));
else
    selectedcell = [goodcell1 goodcell2];
    peak_selected = [peak_selected1;peak_selected2];
    h0 = plotpeakselected_new(peak_selected);
    figname=sprintf('selected_%dpairs',numel(goodcell1));
end

save( [filename(1:end-4) suffix '.mat'],'peak_selected','selectedcell')
save( [figname suffix '.mat'],'peak_selected','selectedcell')
saveas(h0,[ figname suffix '.fig'])
saveas(h0,[ figname suffix '.png'])