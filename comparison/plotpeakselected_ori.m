function h0 = plotpeakselected_ori(peak_selected)
%peak_selected: status(4:day1still,day1run,day2still,day2run)*nSteps*ncell

ncell = size(peak_selected,3);
nSteps = size(peak_selected,2);
nSteps1 = nSteps - mod(nSteps,2);
nOris1 = nSteps1/2;
nOris = nOris1+(nSteps-nSteps1);

[max_amp,pref_ori0]=nanmax(nanmax(peak_selected(:,1:end-1,:)));
pref_ori = pref_ori0;
pref_ori(pref_ori>nOris1) = pref_ori(pref_ori>nOris1)-nOris1;

% ori = squeeze(nanmean(reshape(peak_selected(:,1:nSteps1,:),4,nOris1,2,ncell),3));
% if nOris
%     ori(:,nOris,:) = peak_selected(:,nSteps,:);
% end
% [~,pref_ori] = max(max(ori));
%%
resp1(1,:) = calgOSI(peak_selected(1,1:end-1,:));
resp1(2,:) = calgOSI(peak_selected(2,1:end-1,:));
resp2(1,:) = calgOSI(peak_selected(3,1:end-1,:));
resp2(2,:) = calgOSI(peak_selected(4,1:end-1,:));
resp = cat(1,resp1,resp2);
% pref_ori = nanmax(pref_ori0);
% pref_ori(pref_ori>nOris1) = pref_ori(pref_ori>nOris1)-nOris1;
thr = 0.05;
resp(:,max_amp<thr) = NaN;
resp(:,sum(isnan(resp))>0) = NaN;
resp1 =resp(1:2,:);
resp2 =resp(3:4,:);

% resp2(isinf(resp2)) = NaN;
% resp1(isinf(resp1)) = NaN;
% resp2(resp2>prctile(resp2(:),99)) =NaN;
% resp1(resp1>prctile(resp1(:),99)) =NaN;

for i=1:nOris
    peak1(:,i) = nanmean(resp1(:,pref_ori==i),2);
    peak2(:,i) = nanmean(resp2(:,pref_ori==i),2);
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


figurename=sprintf('%d cells with response greater than %.2f',sum(~sum(isnan(resp))),thr);
h0=figure('Position',[0 0 1200 600],'Name',figurename);
hold on
ymin = min(resp1(:));
ymax = max(resp1(:));
for k=1:nOris
    try
        plot([k k+.4],squeeze(resp1(:,pref_ori==k)),'o-','Color', [.5 .5 .5 .5])
        plot([k+.1 k+.5],squeeze(resp2(:,pref_ori==k)),'x-','Color', [.5 .5 .5 .5])
    end
    plot([k k+.4],squeeze(peak1(:,k)),'ko-','MarkerSize',8,'LineWidth',3)
    plot([k+.1 k+.5],squeeze(peak2(:,k)),'kx--','MarkerSize',8,'LineWidth',3)

end

% title(['Response@' names{i}])
legend('Day0','Day7','Location','northeast')
legend('boxoff')
axis tight
xlim([0.5 nOris+1])
% ylim([0 10]);
ax=gca;
ax.XTick = 1:nOris;
ax.XTickLabel = names;
xlabel('Preferred Orientation')
ylabel('OSI')


