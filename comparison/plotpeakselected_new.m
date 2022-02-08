function h0 = plotpeakselected_new(peak_selected)
%peak_selected: status(4:day1still,day1run,day2still,day2run)*nSteps*ncell

ncell = size(peak_selected,3);
nSteps = size(peak_selected,2);
nSteps1 = nSteps - mod(nSteps,2);
nOris1 = nSteps1/2;
nOris = nOris1+(nSteps-nSteps1);
% nOris1 = nSteps1;
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
%%
for i=1:ncell
    amp(i) = ori(2,pref_ori(i),i);
    resp1(:,:,i) = ori(1:2,:,i)/amp(i);
    resp2(:,:,i) = ori(3:4,:,i)/amp(i);
end
resp1 = real(log2(resp1));
resp2 = real(log2(resp2));

thr = 0.05;
% resp (:,:,normal<thr) = NaN;
resp1(:,:,max_amp<thr) = NaN;
resp2(:,:,max_amp<thr) = NaN;
resp2(isinf(resp2)) = NaN;
resp1(isinf(resp1)) = NaN;
resp1(isnan(resp2)) =NaN;
resp2(isnan(resp1)) =NaN;
thr_up = prctile( [resp1(:);resp2(:)],99);
resp2(resp2>thr_up) =NaN;
resp1(resp1>thr_up) =NaN;


for i=1:nOris
    peak1(:,:,i) = nanmean(resp1(:,:,pref_ori==i),3);
    peak2(:,:,i) = nanmean(resp2(:,:,pref_ori==i),3);
end

%%
h = (1:nOris)/nOris;
s = ones(1,nOris);
v = ones(1,nOris);
clut =hsv2rgb(h,s,v);

for i=1:nOris1
    names{i} = sprintf('%.1f',(i-1)/nSteps1*360);
end
if nOris>nOris1
    names{nOris} = 'blk';
end


figurename=sprintf('%d cells with response greater than %.2f',sum(~sum(isnan(sum(resp1)))),thr);
h0=figure('Position',[0 0 1200 600],'Name',figurename);
for i=1:nOris
    subplot(nOris1/2+1,2,i);hold on
    for k=1:nOris
  plot([k k+.4],squeeze(peak1(:,i,k)),'o-','MarkerSize',8,'LineWidth',3,'Color',clut(1,i,:))
       plot([k+.1 k+.5],squeeze(peak2(:,i,k)),'x--','MarkerSize',8,'LineWidth',3,'Color',clut(1,i,:))
        try
        plot([k k+.4],squeeze(resp1(:,i,pref_ori==k)),'o-','LineWidth',.3,'Color', [.5 .5 .5 .5])
        plot([k+.1 k+.5],squeeze(resp2(:,i,pref_ori==k)),'x-','LineWidth',.3,'Color', [.5 .5 .5 .5])
       end
       plot([k k+.4],squeeze(peak1(:,i,k)),'o-','MarkerSize',8,'LineWidth',3,'Color',clut(1,i,:))
       plot([k+.1 k+.5],squeeze(peak2(:,i,k)),'x--','MarkerSize',8,'LineWidth',3,'Color',clut(1,i,:))
 
    end
    
    title(['Response@' names{i}])
    legend('Day0','Day7','Location','northeast')
    legend('boxoff')
    axis tight
    xlim([0.5 nOris+1])
    ylim([.33 thr_up]);
    ax=gca;
%    ax.YScale = 'linear';
    ax.YScale = 'log';
    ax.XTick = 1:nOris;
    ax.XTickLabel = names;
end

subplot(nOris1/2+1,2,nOris1+2);hold on
% figure;hold on
for k=1:nOris
    plot([k k+.4],squeeze(peak1(:,k,k)),'o-','MarkerSize',8,'LineWidth',3,'Color',clut(1,k,:))
    plot([k+.1 k+.5],squeeze(peak2(:,k,k)),'x--','MarkerSize',8,'LineWidth',3,'Color',clut(1,k,:))
end
ymax = max([peak1(:);peak2(:)]);
xlim([0.5 nOris+1])
ylim( [0 ceil(ymax)])
xlabel('Preferred Orientation')
ylabel('Amplitude Ratio')
ax=gca;
ax.XTick = 1:nOris;
ax.XTickLabel = names;

