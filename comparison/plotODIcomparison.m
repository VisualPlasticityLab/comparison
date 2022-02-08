%function plotODIcomparison
clear all;clc
[day1peakf,path1]=uigetfile('*all.mat','select ALL cell file for day1');
[day1pairf,path1a]=uigetfile([path1,'*as+ms*.mat'],'select picked pairs for day1');
[day2peakf,path2]=uigetfile('*all.mat','select ALL cell file for day2');
[day2pairf,path2a]=uigetfile([path2 '*as+ms*.mat'],'select picked pairs for day2');
[prefilename,path]=uigetfile('*.mat','Pick the pre-matched pairs?');
% ODI0(jj) = odi;
% ODI1(:,jj) = [odi_s  odi_r]; %using same pref_ori
% PEAK0(jj) =a;
% PEAK00(:,jj)= [eye1.peak(:,pref_ori,nth1) eye1.peakS(:,pref_ori,nth1) eye1.peakR(:,pref_ori,nth1)...
%        eye2.peak(:,pref_ori,nth2) eye2.peakS(:,pref_ori,nth2) eye2.peakR(:,pref_ori,nth2)]; % use combined pref_ori during after/before 
day1peak=load(fullfile(path1,day1peakf),'PEAK0','PEAK00','ODI0','ODI1');
day2peak=load(fullfile(path2,day2peakf),'PEAK0','PEAK00','ODI0','ODI1');
slt1=load(fullfile(path1a,day1pairf));
slt2=load(fullfile(path2a,day2pairf));
pre = load(fullfile(fullfile(path,prefilename)),'pair1','pair2');
%% PICK cells at least have a good response in one day
[~,IP1] = intersect(pre.pair1,find(slt1.pair0&slt1.pair1&slt1.pair2))
[~,IP2] = intersect(pre.pair2,find(slt2.pair0&slt2.pair1&slt2.pair2))
eitherday = union(IP1,IP2);
numcell = num2str(numel(eitherday));

pairday1 = pre.pair1(eitherday);
pairday2 = pre.pair2(eitherday);
% level = 6;
% peak = discretize(day1peak.PEAK0(pair1),level)/level./day1peak.PEAK0(pair1);

% STILL& RUN condition
clear ODI1s ODI1r
ODI1s(1,:) = day1peak.ODI1(1,pairday1);
ODI1s(2,:) = day2peak.ODI1(1,pairday2);

ODI1r(1,:) = day1peak.ODI1(2,pairday1);
ODI1r(2,:) = day2peak.ODI1(2,pairday2);

pos_ = strfind(prefilename,'_');
filenm = [prefilename(1:pos_(end)) 'selected' numcell];

%% plots all overlapped
load([path prefilename(1:pos_(end)) 'aligned.mat'])
fig0 = figure('Position',[0 0 2000 1200]);
SIZ = size(fig2img);
imshowpair(fig1imgnew,fig2img)
hold on;
for i=pre.pair1
    contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,i),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[1 0 0])
end
%2nd image, green
for i=pre.pair2
    contour(1:SIZ(2),1:SIZ(1),reshape(F2Contour(:,i),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[0 1 0])
end
%%
saveas(fig0,[filenm '.fig'])
%%
sets = 2; %define how many columns
h(2)=figure('Position',[200 200 300*sets 500]);

subplot(2,sets,1);title('STILL','FontSize',12)
hold on
plot([-1 1],[-1 1],'--','Color',[0.5 0.5 .5 .5])
% scatter(ODI1s(1,:),ODI1s(2,:),peak*50,'ko','filled'); 
scatter(ODI1s(1,:),ODI1s(2,:),'ko','filled'); 
scatter(nanmean(ODI1s(1,:)),nanmean(ODI1s(2,:)),100,'kx');
alpha(.3)
xlabel(['ODIbefore \mu= ' sprintf('%.2f',nanmean(ODI1s(1,:)))])
ylabel(['ODIafter \mu= ' sprintf('%.2f',nanmean(ODI1s(2,:)))])
axis square

subplot(2,sets,sets+1);hold on;
histplt2(ODI1s(1,:),ODI1s(2,:),-1.0:2/11:1.0);
legend('before','after','Location','Northwest')
legend('Boxoff')
xlabel('ODI@best direction')
ylabel('Cell count')
title(['N=' numcell],'FontSize',12)
axis square

% RUN  condition
subplot(2,sets,2);title('RUN','FontSize',12)
hold on
plot([-1 1],[-1 1],'--','Color',[0.5 0.5 .5 .5])
% scatter(ODI1r(1,:),ODI1r(2,:),peak*50,'ko','filled'); 
scatter(ODI1r(1,:),ODI1r(2,:),'ko','filled'); 
scatter(nanmean(ODI1r(1,:)),nanmean(ODI1r(2,:)),100,'kx');
alpha(.3)

xlabel(['ODIbefore \mu= ' sprintf('%.2f',nanmean(ODI1r(1,:)))])
ylabel(['ODIafter \mu= ' sprintf('%.2f',nanmean(ODI1r(2,:)))])
axis square

subplot(2,sets,sets+2);hold on;
histplt2(ODI1r(1,:),ODI1r(2,:),-1.0:2/11:1.0);
legend('before','after','Location','Northwest')
legend('Boxoff')
xlabel('ODI@best direction')
ylabel('Cell count')
title(['N=' numcell],'FontSize',12)
axis square

%%
saveas(h(2),[filenm '.png'])
saveas(h(2),[filenm '.fig'])

save(filenm,'ODI1r','ODI1s','IP1','IP2','pairday1','pairday2')

%%
