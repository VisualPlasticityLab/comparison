function h0 = plotpeakselected_drn(peak_selected)
%peak_selected: status(4:day1still,day1run,day2still,day2run)*nSteps*ncell

ncell = size(peak_selected,3);
nSteps = size(peak_selected,2);
nSteps1 = nSteps - mod(nSteps,2);
nDrns = nSteps1+(nSteps-nSteps1);

% ori = (peak_selected(:,1:nSteps1/2,:)+peak_selected(:,nSteps1/2+1:nSteps1,:))/2;
drn = peak_selected;
[~,pref_ori] = max(max(drn));
%%
for i=1:ncell
    amp(i) = drn(1,pref_ori(i),i);
    resp1(:,:,i) = drn(1:2,:,i)/amp(i);
    resp2(:,:,i) = drn(3:4,:,i)/amp(i);
end
thr = 0.05;
% resp (:,:,normal<thr) = NaN;

for i=1:nDrns
    peak1(:,:,i) = nanmedian(resp1(:,:,pref_ori==i),3);
    peak2(:,:,i) = nanmedian(resp2(:,:,pref_ori==i),3);
end

%%


h = (1:nDrns)/nDrns;
s = ones(1,nDrns);
v = ones(1,nDrns);
clut =hsv2rgb(h,s,v);

for i=1:nSteps1
    names{i} = sprintf('%.1f',(i-1)/nSteps1*360);
end
if nDrns>nSteps1
    names{nDrns} = 'blk';
end


figurename=sprintf('%d cells with response greater than %.2f',sum(amp>thr),thr);
h0=figure('Position',[0 0 1200 1200],'Name',figurename);
ymin = min(resp1(:));
ymax = max(resp1(:));
for i=1:nDrns
    subplot(nSteps1/2+1,2,i);hold on
    for k=1:nDrns
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
    xlim([0.5 nDrns+1])
    ax=gca;
    ax.XTick = 1:nDrns;
    ax.XTickLabel = names;
end

subplot(nSteps1/2+1,2,nDrns+1);hold on
for k=1:nDrns
    plot([k k+.4],squeeze(peak1(:,k,k)),'o-','MarkerSize',8,'LineWidth',3,'Color',clut(1,k,:))
    plot([k+.1 k+.5],squeeze(peak2(:,k,k)),'x--','MarkerSize',8,'LineWidth',3,'Color',clut(1,k,:))
end
xlim([0.5 nDrns+1])
xlabel('Preferred Direction')
ax=gca;
ax.XTick = 1:nDrns;
ax.XTickLabel = names;

