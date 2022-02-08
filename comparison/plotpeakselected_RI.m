function h0 = plotpeakselected_RI(peak_selected)
%peqk_selected: status(4:day1still,day1run,day2still,day2run)*nSteps*ncell

ncell = size(peak_selected,3);
nSteps = size(peak_selected,2);
nSteps1 = nSteps - mod(nSteps,2);
nOris1 = nSteps1/2;
nOris = nOris1+(nSteps-nSteps1);
% nOris1 = nSteps1;edi
% nOris = nSteps;

[max_amp,pref_ori0]=nanmax(nanmax(peak_selected));
pref_ori = pref_ori0;
pref_ori(pref_ori>nOris1) = pref_ori(pref_ori>nOris1)-nOris1;

ori = squeeze(nanmean(reshape(peak_selected(:,1:nSteps1,:),4,nOris1,2,ncell),3));
% ori = peak_selected;
if nOris
    ori(:,nOris,:) = peak_selected(:,nSteps,:);
end
% [~,pref_ori] = max(max(ori));

for i=1:ncell
    resp(1,:,i) = real(log2(ori(2,pref_ori(i),i)/ori(1,pref_ori(i),i)));
    resp(2,:,i) = real(log2(ori(4,pref_ori(i),i)/ori(3,pref_ori(i),i)));
end
%% Only max_amp > thr and not isinf/isnan in any days
thr = 0.05;
% resp (:,:,normal<thr) = NaN;
resp(:,:,max_amp<thr) = NaN;
resp(isinf(resp)) = NaN;
% thr_up = prctile( resp(:),[1 99]);
% resp(resp>thr_up(2)) =NaN;
% resp(resp<thr_up(1)) =NaN;
% resp(:,:,isnan(sum(resp,1))) = NaN;

resp1=resp(1,:,:) ;
resp2=resp(2,:,:) ;


for i=1:nOris
    selected=(pref_ori==i)&(~isnan(resp1)&~isnan(resp2));
    count(i) = sum(selected);
    peak1(:,:,i) = mean(resp1(:,:,selected),3);
    peak2(:,:,i) = mean(resp2(:,:,selected),3);
    var1(:,:,i) = std(squeeze(resp1(:,:,selected)))/sqrt(count(i));
    var2(:,:,i) = std(squeeze(resp2(:,:,selected)))/sqrt(count(i));
end

%%
h = (1:nOris)/nOris;
s = ones(1,nOris);
v = ones(1,nOris);
clut =hsv2rgb(h,s,v);

for i=1:nOris1
    names{i} = sprintf('%d/%d',(i-1)/nSteps1*360,(i-1)/nSteps1*360+180);
end
if nOris>nOris1
    names{nOris} = 'blk';
end

%% plot running index
figurename=sprintf('Modulation Index :%d cells with response greater than %.2f',sum(count),thr);
h0=figure('Position',[200 200 1200 300],'Name',figurename);hold on;
i=1;
% subplot(nOris1/2+1,2,i);hold on
for k=1:nOris
    errorbar([k;k+.4],squeeze(cat(1,peak1(:,i,k),peak2(:,i,k))),squeeze(cat(1,var1(:,i,k),var2(:,i,k))),'o-','MarkerSize',4,'LineWidth',2,'Color',clut(1,i+2,:))
    %        plot([k+.1 k+.5],squeeze(peak2(:,i,k)),'x--','MarkerSize',8,'LineWidth',3,'Color',clut(1,i,:))
    try
        plot([k;k+.4],squeeze(cat(1,resp1(:,i,pref_ori==k),resp2(:,i,pref_ori==k))),'o-','LineWidth',.3,'Color', [.5 .5 .5 .5])
    end
        errorbar([k;k+.4],squeeze(cat(1,peak1(:,i,k),peak2(:,i,k))),squeeze(cat(1,var1(:,i,k),var2(:,i,k))),'o-','MarkerSize',4,'LineWidth',2,'Color',clut(1,i+2,:))

    [T(k),P(k)]=ttest(resp1(:,i,pref_ori==k),resp2(:,i,pref_ori==k));
end
legend('Day0--Day7','Location','northeast')
legend('boxoff')
axis tight

ax=gca;

ax.YScale = 'linear';
ylim([-5 5]);
% 
% ax.YScale = 'log';
% ylim([.025 40]);
% ax.YTick = 10.^(-1.6:1.6:1.6);
% ax.YTickLabel = [.025,1, 40];

xlim([0.5 nOris+1])

ax.XTick = (1:nOris)+.25;
ax.XTickLabel = names;

xlabel('Preferred Orientation');
ylabel('Ratio(run/still)');
title(figurename)

[T,P]



%
% subplot(nOris1/2+1,2,nOris1+2);hold on
% % figure;hold on
% for k=1:nOris
%     plot([k k+.4],squeeze(peak1(:,k,k)),'o-','MarkerSize',8,'LineWidth',3,'Color',clut(1,k,:))
%     plot([k+.1 k+.5],squeeze(peak2(:,k,k)),'x--','MarkerSize',8,'LineWidth',3,'Color',clut(1,k,:))
% end
% ymax = max([peak1(:);peak2(:)]);
% xlim([0.5 nOris+1])
% ylim( [0 ceil(ymax)])
% xlabel('Preferred Orientation')
% ylabel('Amplitude Ratio')
% ax=gca;
% ax.XTick = 1:nOris;
% ax.XTickLabel = names;

%%
% figurename=sprintf('%d cells with response greater than %.2f',sum(~sum(isnan(sum(resp1)))),thr);
% h0=figure('Position',[0 0 1200 600],'Name',figurename);
% for i=1:nOris
%     subplot(nOris1/2+1,2,i);hold on
%     for k=1:nOris
%        plot([k k+.4],squeeze(peak1(:,i,k)),'o-','MarkerSize',8,'LineWidth',3,'Color',clut(1,i,:))
%        plot([k+.1 k+.5],squeeze(peak2(:,i,k)),'x--','MarkerSize',8,'LineWidth',3,'Color',clut(1,i,:))
%        try
%         plot([k k+.4],squeeze(resp1(:,i,pref_ori==k)),'o-','LineWidth',.3,'Color', [.5 .5 .5 .5])
%         plot([k+.1 k+.5],squeeze(resp2(:,i,pref_ori==k)),'x-','LineWidth',.3,'Color', [.5 .5 .5 .5])
%        end
%        plot([k k+.4],squeeze(peak1(:,i,k)),'o-','MarkerSize',8,'LineWidth',3,'Color',clut(1,i,:))
%        plot([k+.1 k+.5],squeeze(peak2(:,i,k)),'x--','MarkerSize',8,'LineWidth',3,'Color',clut(1,i,:))
%
%     end
%
%     title(['Response@' names{i}])
%     legend('Day0','Day7','Location','northeast')
%     legend('boxoff')
%     axis tight
%     xlim([0.5 nOris+1])
%     ylim([.33 thr_up]);
%     ax=gca;
% %    ax.YScale = 'linear';
%     ax.YScale = 'log';
%     ax.XTick = 1:nOris;
%     ax.XTickLabel = names;
% end
%
% subplot(nOris1/2+1,2,nOris1+2);hold on
% % figure;hold on
% for k=1:nOris
%     plot([k k+.4],squeeze(peak1(:,k,k)),'o-','MarkerSize',8,'LineWidth',3,'Color',clut(1,k,:))
%     plot([k+.1 k+.5],squeeze(peak2(:,k,k)),'x--','MarkerSize',8,'LineWidth',3,'Color',clut(1,k,:))
% end
% ymax = max([peak1(:);peak2(:)]);
% xlim([0.5 nOris+1])
% ylim( [0 ceil(ymax)])
% xlabel('Preferred Orientation')
% ylabel('Amplitude Ratio')
% ax=gca;
% ax.XTick = 1:nOris;
% ax.XTickLabel = names;


% %%
% figurename=sprintf('%d cells with response greater than %.2f',sum(~sum(isnan(sum(resp1)))),thr);
% h0=figure('Position',[0 0 1200 600],'Name',figurename);
% for i=1:nOris
%     subplot(nOris1/2+1,2,i);hold on
%     for k=1:nOris
%        plot([k k+.4],squeeze(peak1(:,i,k)),'o-','MarkerSize',8,'LineWidth',3,'Color',clut(1,i,:))
%        plot([k+.1 k+.5],squeeze(peak2(:,i,k)),'x--','MarkerSize',8,'LineWidth',3,'Color',clut(1,i,:))
%        try
%         plot([k k+.4],squeeze(resp1(:,i,pref_ori==k)),'o-','LineWidth',.3,'Color', [.5 .5 .5 .5])
%         plot([k+.1 k+.5],squeeze(resp2(:,i,pref_ori==k)),'x-','LineWidth',.3,'Color', [.5 .5 .5 .5])
%        end
%        plot([k k+.4],squeeze(peak1(:,i,k)),'o-','MarkerSize',8,'LineWidth',3,'Color',clut(1,i,:))
%        plot([k+.1 k+.5],squeeze(peak2(:,i,k)),'x--','MarkerSize',8,'LineWidth',3,'Color',clut(1,i,:))
%
%     end
%
%     title(['Response@' names{i}])
%     legend('Day0','Day7','Location','northeast')
%     legend('boxoff')
%     axis tight
%     xlim([0.5 nOris+1])
%     ylim([.33 thr_up]);
%     ax=gca;
% %    ax.YScale = 'linear';
%     ax.YScale = 'log';
%     ax.XTick = 1:nOris;
%     ax.XTickLabel = names;
% end
%
% subplot(nOris1/2+1,2,nOris1+2);hold on
% % figure;hold on
% for k=1:nOris
%     plot([k k+.4],squeeze(peak1(:,k,k)),'o-','MarkerSize',8,'LineWidth',3,'Color',clut(1,k,:))
%     plot([k+.1 k+.5],squeeze(peak2(:,k,k)),'x--','MarkerSize',8,'LineWidth',3,'Color',clut(1,k,:))
% end
% ymax = max([peak1(:);peak2(:)]);
% xlim([0.5 nOris+1])
% ylim( [0 ceil(ymax)])
% xlabel('Preferred Orientation')
% ylabel('Amplitude Ratio')
% ax=gca;
% ax.XTick = 1:nOris;
% ax.XTickLabel = names;
%
