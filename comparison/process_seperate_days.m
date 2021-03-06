%% Load Day1 and Day2 traces and peakSI info
% [day1sigf,path]=uigetfile('*signals.mat','select trace for day1eye');
% day1sig=load(fullfile(path,day1sigf))
% day1.sigF= permute(day1sig.sigF,[4 3 2 1]);
[day1peakf,path]=uigetfile('peakSI.mat','select peakSI matfile for day1eye');
day1peak=load(fullfile(path,day1peakf));
day1.finalvalue = permute(day1peak.SI.peakS,[3 2 1]);
day1.stdofeachvalue = permute(day1peak.SI.errorS,[3 2 1]);
day1.finalvalue2 = permute(day1peak.SI.peakR,[3 2 1]);
day1.stdofeachvalue2 = permute(day1peak.SI.errorR,[3 2 1]);
day1.sigF= permute(day1peak.sigF,[4 3 2 1]);
% [day2sigf,path]=uigetfile('*signals.mat','open trace for day2eye');
% day2sig=load(fullfile(path,day2sigf));
% day2.sigF= permute(day2sig.sigF,[4 3 2 1]);
[day2peakf,path]=uigetfile('peakSI.mat','select peakSI matfile for day2eye');
day2peak=load(fullfile(path,day2peakf));
day2.finalvalue = permute(day2peak.SI.peakS,[3 2 1]);
day2.stdofeachvalue = permute(day2peak.SI.errorS,[3 2 1]);
day2.finalvalue2 = permute(day2peak.SI.peakR,[3 2 1]);
day2.stdofeachvalue2 = permute(day2peak.SI.errorR,[3 2 1]);
day2.sigF= permute(day2peak.sigF,[4 3 2 1]);
%% sigF:ncells,nSteps, numberofcycles,framestocapture

nSteps=size(day1.sigF,2);
framestocapture=min (size(day1.sigF,4), size(day2.sigF,4));
nSteps1=nSteps-mod(nSteps,2);

goodcell1 = 1:size(day1.finalvalue,1);
goodcell2 = 1:size(day2.finalvalue,1);
%%
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

for i=1:size(peak_selected1,3)
%     normal(i)=max(max(peak_selected(:,:,i)));
     normal(i) = max(max(peak_selected1(:,:,i)));
    resp1(1:2,:,i) = peak_selected1(1:2,:,i)/normal(i);
end

for i=1:size(peak_selected2,3)
%     normal(i)=max(max(peak_selected(:,:,i)));
    normal2(i) = max(max(peak_selected2(:,:,i)));
    resp2(1:2,:,i) = peak_selected2(1:2,:,i)/normal2(i);
end
%%
thr = 0;
% resp (:,:,normal<thr) = NaN;
% resp (:,:,8) = NaN;
% resp (:,:,12) = NaN;

figurename=sprintf('%d cells with response greater than %.2f',sum(normal>thr),thr);
figure('Position',[0 0 1200 1200],'Name',figurename);
ymin = min(resp1(:));
ymax = max(resp1(:));

traces0= squeeze(resp1(:,1,:));
traces45= squeeze(resp1(:,6,:));
traces90=squeeze(resp1(:,3,:));
traces135 = squeeze(resp1(:,4,:));

traces0R= squeeze(resp2(:,1,:));
traces45R= squeeze(resp2(:,6,:));
traces90R= squeeze(resp2(:,3,:));
traces135R= squeeze(resp2(:,4,:));

traces90 = resp1(:,7,:)

resp2(:,7,:)
resp2(:,8,:)
resp1(:,8,:)

subplot(3,2,1);hold on

plot(1:2,traces0,'x-','Color', [.5 .5 .5])
plot(1:2,nanmean(traces0,2),'o-','MarkerSize',8,'LineWidth',3)
plot(3:4,traces0R,'x-','Color', [.5 .5 .5])
plot(3:4,nanmean(traces0R,2),'o-','MarkerSize',8,'LineWidth',3)

title('Response ORI=0')
ylim([ymin ymax])

subplot(3,2,2);hold on
% traces45= squeeze((resp1(:,2,:)+resp1(:,6,:))/2);
plot(1:2,traces45,'x-','Color', [.5 .5 .5])
plot(1:2,nanmean(traces45,2),'o-','MarkerSize',8,'LineWidth',3)
% traces45R= squeeze((resp2(:,2,:)+resp2(:,6,:))/2);
plot(3:4,traces45R,'x-','Color', [.5 .5 .5])
plot(3:4,nanmean(traces45R,2),'o-','MarkerSize',8,'LineWidth',3)

title('Response ORI=\pi/4')
ylim([ymin ymax])

subplot(3,2,3);hold on
plot(1:2,traces90,'x-','Color', [.5 .5 .5])
plot(1:2,nanmean(traces90,2),'o-','MarkerSize',8,'LineWidth',3)
plot(3:4,traces90R,'x-','Color', [.5 .5 .5])
plot(3:4,nanmean(traces90R,2),'o-','MarkerSize',8,'LineWidth',3)
title('Response ORI=\pi/2')
ylim([ymin ymax])

subplot(3,2,4);hold on

plot(1:2,traces135,'x-','Color', [.5 .5 .5])
plot(1:2,nanmean(traces135,2),'o-','MarkerSize',8,'LineWidth',3)
plot(3:4,traces135R,'x-','Color', [.5 .5 .5])
plot(3:4,nanmean(traces135R,2),'o-','MarkerSize',8,'LineWidth',3)
title('Response ORI=3/4\pi')
ylim([ymin ymax])

subplot(3,2,5);hold on
tracesblk=squeeze(resp1(:,9,:));
plot(1:2,tracesblk,'x-','Color', [.5 .5 .5])
plot(1:2,nanmean(tracesblk,2),'o-','MarkerSize',8,'LineWidth',3)
tracesblkR= squeeze(resp2(:,9,:));
plot(3:4,tracesblkR,'x-','Color', [.5 .5 .5])
plot(3:4,nanmean(tracesblkR,2),'o-','MarkerSize',8,'LineWidth',3)

title('Response blank')
ylim([ymin ymax])

subplot(3,2,6);hold on
plot(1:2,nanmean(traces0,2),'o-','MarkerSize',8,'LineWidth',3)
plot(1:2,nanmean(traces45,2),'o-','MarkerSize',8,'LineWidth',3)
plot(1:2,nanmean(traces90,2),'o-','MarkerSize',8,'LineWidth',3)
plot(1:2,nanmean(traces135,2),'o-','MarkerSize',8,'LineWidth',3)
plot(1:2,nanmean(tracesblk,2),'o-','MarkerSize',8,'LineWidth',3)

plot(3:4,nanmean(traces0R,2),'o-','MarkerSize',8,'LineWidth',3)
plot(3:4,nanmean(traces45R,2),'o-','MarkerSize',8,'LineWidth',3)
plot(3:4,nanmean(traces90R,2),'o-','MarkerSize',8,'LineWidth',3)
plot(3:4,nanmean(traces135R,2),'o-','MarkerSize',8,'LineWidth',3)
plot(3:4,nanmean(tracesblkR,2),'o-','MarkerSize',8,'LineWidth',3)

legend('0','\pi/4','\pi/2','3/4\pi','blank','0R','\pi/4R','\pi/2R','3/4\piR','blankR''Location','northwest')
legend('boxoff')

% saveas(gcf,[filename '.fig'])
