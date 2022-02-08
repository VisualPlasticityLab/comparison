function h0 = plotpeakselected_direction(peak_selected)
%peqk_selected: status(4:day1still,day1run,day2still,day2run)*nSteps*ncell

ncell = size(peak_selected,3);
nSteps = size(peak_selected,2);
nSteps1 = nSteps - mod(nSteps,2);
nOris1 = nSteps1/2;
nOris = nOris1+(nSteps-nSteps1);

% ori = (peak_selected(:,1:nSteps1/2,:)+peak_selected(:,nSteps1/2+1:nSteps1,:))/2;


ori = nanmean(reshape(peak_selected(:,1:nSteps1,:),4,nOris1,2,ncell),3);
if nOris
    ori(:,nOris,:) = peak_selected(:,nSteps,:);
end
%[amp,pref_ori] = max(max(ori));
[amp,pref_ori] = max(ori(1,:,:));
%%
for i=1:ncell
    resp1(:,:,i) = ori(1:2,:,i)/amp(i);
    resp2(:,:,i) = ori(3:4,:,i)/amp(i);
end
thr = 0;
% resp (:,:,normal<thr) = NaN;

for i=1:nOris
    peak1(:,:,i) = nanmedian(resp1(:,:,pref_ori==i),3);
    peak2(:,:,i) = nanmedian(resp2(:,:,pref_ori==i),3);
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


figurename=sprintf('%d cells with response greater than %.2f',sum(amp>thr),thr);
h0=figure('Position',[0 0 1200 1200],'Name',figurename);
ymin = min(resp1(:));
ymax = max(resp1(:));
for i=1:nOris
    subplot(nSteps1/4+1,2,i);hold on
    for k=1:nOris
        plot([k k+.4],squeeze(peak1(:,i,k)),'o-','MarkerSize',8,'LineWidth',3,'Color',clut(1,i,:))
        plot([k+.1 k+.5],squeeze(peak2(:,i,k)),'x--','MarkerSize',8,'LineWidth',3,'Color',clut(1,i,:))
        try
        plot([k k+.4],squeeze(resp1(:,i,pref_ori==k)),'o-','Color', [.5 .5 .5 .5])
        plot([k+.1 k+.5],squeeze(resp2(:,i,pref_ori==k)),'x-','Color', [.5 .5 .5 .5])
        end
    end
    
    title(['Response@' names{i}])
    legend('Day0','Day7','Location','northeast')
    legend('boxoff')
    axis tight
    xlim([0.5 nOris+1])
    ax=gca;
    ax.XTick = 1:nOris;
    ax.XTickLabel = names;
end

subplot(nSteps1/4+1,2,nOris+1);hold on
for k=1:nOris
    plot([k k+.4],squeeze(peak1(:,k,k)),'o-','MarkerSize',8,'LineWidth',3,'Color',clut(1,k,:))
    plot([k+.1 k+.5],squeeze(peak2(:,k,k)),'x--','MarkerSize',8,'LineWidth',3,'Color',clut(1,k,:))
end
xlim([0.5 nOris+1])
xlabel('Preferred Orientation')
ax=gca;
ax.XTick = 1:nOris;
ax.XTickLabel = names;

