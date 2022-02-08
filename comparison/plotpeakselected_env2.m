function h0 = plotpeakselected_env2(peak_selected)
%% similiar to plotpeakselected_env except calculating direction
%% peak_selected: status(2or4:day1still,day1run,day2still,day2run)*nSteps*ncell
thr = 0.05;
N = [3,4];
Nstats = 'Day1';

peak_selected(peak_selected<0) = 0;
peak_selected(:,:,max(max(peak_selected))<thr) = [];

nStatus = size(peak_selected,1);
nSteps = size(peak_selected,2);
nSteps1 = nSteps - mod(nSteps,2);
ori = nSteps1/2;

     [amp_max1,pref_dir1]=nanmax(nanmax(peak_selected(3:4,:,:)));
     [amp_max0,pref_dir0]=nanmax(nanmax(peak_selected(1:2,:,:)));
     [amp_min1,~]=nanmin(nanmin(peak_selected(3:4,:,:)));
     [amp_min0,~]=nanmin(nanmin(peak_selected(1:2,:,:)));


 % use Day0 or Day1 as the preferred direntation
 if numel(N)>1
     [amp_max,pref_dir]=max(nanmax(peak_selected(N,:,:)));
     [amp_min,~]=min(nanmin(peak_selected(N,:,:)));

 else
     
     [amp_max,pref_dir]=nanmax(peak_selected(N,:,:));
     [amp_min,~]=nanmin(peak_selected(N,:,:));
    
 end
% th_max=prctile(amp_max(:),[2 98]);
% resp(:,:,amp_max<th_max(1)|amp_max>th_max(2)) = [];
% 
% th_min=prctile(amp_min,[2 98]);
% amp_min(amp_min<th_min(1)|amp_min>th_min(2))= NaN;

%% Only max_amp > thr and not isinf/isnan in any days
for i=1:size(peak_selected,3)
    resp(:,:,i) = (peak_selected(:,:,i)-amp_min(i))/(amp_max(i)-amp_min(i));
%     resp(1:2,:,i) = (peak_selected(1:2,:,i)-amp_min0(i))/(amp_max0(i)-amp_min0(i));
%     resp(3:4,:,i) = (peak_selected(3:4,:,i)-amp_min1(i))/(amp_max1(i)-amp_min1(i));

end
% resp = real(log2(resp+1));

% resp (:,:,normal<thr) = NaN;
resp(isinf(resp)) = NaN;
resp(resp<0) == 0 ;
thr_up = prctile(resp(:),[1 99]);
resp(resp<thr_up(1)|resp>thr_up(2)) =NaN;
%resp size: 4*nSteps*cells
%peak size: 4*nSteps*pref_dir
%peak_dir size: 4*pref_dir

%%  Set up plotting labels and color schemes
% colortype = nOris;

figurename=sprintf('Cells%d Amp%.2f+ Normed to %s',size(peak_selected,3),thr,Nstats);
colortype = 2;
h = (1:colortype)'/colortype;
s = ones(colortype,1);
v = ones(colortype,1);
clut =squeeze(hsv2rgb(h,s,v));
clut(3:4,:) =hsv2rgb(h,s,v)/2;

[names,names_dir] = num2ori(nSteps);

%resp: (day1still,day1run,day2still,day2run)*nSteps*ncell
for i=1:nSteps
    selected{i}=find(pref_dir==i);%&(~isnan(sum(sum(resp))));
    count(i) = numel(selected{i});
    pref_resp{i} = squeeze(resp(:,i,selected{i}));
    peak(:,:,i) = nanmean(resp(:,:,selected{i}),3);
    var(:,:,i) = nanstd(resp(:,:,selected{i}),0,3)/sqrt(count(i));
    peak_po(:,i)= nanmean(pref_resp{i},2);
    var_po(:,i) = nanstd(pref_resp{i},1,2)/sqrt(count(i));
end
%% plot Overall changes
figure('Name','Population Histogram','Position',[100 800 1800 600]);
subplot(3,3,1);hold on%plot(squeeze(pref_dir0),squeeze(pref_dir1),'o');
% histplt2(squeeze(pref_dir0),squeeze(pref_dir1));
curveplt2(squeeze(pref_dir0),squeeze(pref_dir1),ori);
ax=gca;
ax.XTick = 1:nSteps;
ax.XTickLabel = names_dir;
ylabel('Cell count')
title('preferred direntation');
legend('Day0','Day1','Location','Best');
legend('boxoff')
% subplot(2,2,2);histogram(pref_dir1);title('preferred direction')
% subplot(3,3,2);histplt2(amp_max0,amp_max1,[0:max(amp_max1)/10:max(amp_max1)]);
% title('amp\_max');legend('Day0','Day1');legend('boxoff')
% 
% subplot(3,3,3);%histplt2(amp_min0,amp_min1,[0:max(amp_min0)/10:max(amp_min0)]);
% title('amp\_min');legend('Day0','Day1');legend('boxoff')

% subplot(2,1,2);hold on
% b=bar(1:nSteps,peak_dir');
for k=1:nSteps
    ax(k) = subplot(3,nSteps1+1,k+nSteps1+1);hold on
    %     histplt2(resp(1,k,selected{k}),resp(2,k,selected{k}),linspace(0,1,10));
    %'PlotStyle', 'compact',,'Boxstyle','filled',
    %
%     plot([.15;.3], log2(squeeze(resp([3 4],k,selected{k})./resp(N,k,selected{k}))),'o-','MarkerSize',4,'Color',clut(4,:),'LineWidth',1);
%     plot([-0.15;0.15], log2(squeeze(resp([2 3],k,selected{k})./resp(N,k,selected{k}))),'o-','MarkerSize',4,'Color',clut(3,:),'LineWidth',1);
%     plot([-0.3;-0.15], log2(squeeze(resp([1 2],k,selected{k})./resp(N,k,selected{k}))),'o-','MarkerSize',4,'Color',clut(2,:),'LineWidth',1);
    if count(k) >0
    plot([.15;.3], squeeze(resp([3 4],k,selected{k})),'o-','MarkerSize',4,'Color',clut(4,:),'LineWidth',1);
    plot([-0.15;0.15], squeeze(resp([2 3],k,selected{k})),'o-','MarkerSize',4,'Color',clut(3,:),'LineWidth',1);
    plot([-0.3;-0.15], squeeze(resp([1 2],k,selected{k})),'o-','MarkerSize',4,'Color',clut(2,:),'LineWidth',1);
    if count(k)>1
        boxplot(ax(k),squeeze(resp(:,k,selected{k}))',...
            'Colors',clut,'PlotStyle', 'compact','Boxstyle','outline','Positions',[-.3 -.15 .15 .3],'Labels',{'','Day0','Day1',''});
%         boxplot(ax(k),log2(squeeze(resp(:,k,selected{k})./resp(N,k,selected{k}))'),...
%             'Colors',clut,'PlotStyle', 'compact','Boxstyle','outline','Positions',[-.3 -.15 .15 .3],'Labels',{'','Day0','Day1',''});
    else
        set(gca,'XTick',[-.3 -.15 .15 .3])
        set(gca,'XTickLabel',{'','Day0','Day1',''})
        set(gca,'XTickLabelRotation',90)
    end
    end
    axis tight
    xlim([-.5 .5])
    set(ax,'YTick',[-1 0 1])
    set(ax,'YTickLabel',{'50%', '100%','200%'})
    %xlim([.25 nSteps+1])
    % ylim([-10 10])
    %set(gca,'XTick',[])
    %set(gca,'XTickLabel',{'','Day0','Day1',''})
    title(sprintf('Best dir=%s',names_dir{k}))
end
linkaxes(ax)

for k=1:nSteps
    ax2(k) = subplot(3,nSteps1+1,k+2*nSteps1+2);

%    histogram(pref_dir1(pref_dir0==k),.5:nSteps+.5)
try    histplt2(pref_dir0(pref_dir==k),pref_dir1(pref_dir==k));
    ax2(k).XTick = 1:nSteps;
    ax2(k).XTickLabel = names_dir;
    ylabel('Cell count')
%         set(gca,'XTickLabelRotation',90)
%     axis tight
%     title(sprintf('Best dir=%s',names_dir{k}))

end
end
linkaxes(ax2)

% title('Response peak@preferred direntation')
saveas(gcf,[figurename '_stats.fig'])
saveas(gcf,[figurename '_stats.png'])

%% plot Tuning curve
%resp size: 4*nSteps*cells
%peak size: 4*nSteps*pref_dir
h0=figure('Position',[100 100 1800 600],'Name',figurename);
xaxis = [1:nSteps];% nSteps1+2.5];
fitx = linspace(1,nSteps1,nSteps1*3);


for k=1:nSteps
    if k<nSteps1+1
        ax3(k) = subplot(2,nSteps1/2+1,k+(k>nSteps1/2));
    else
        ax3(k) = subplot(2,nSteps1/2+1,nSteps1/2+1);
    end
    % plot tuning curve @ preferred direntation
    hold on
    for i = 1:nStatus
        errorbar(xaxis,peak(i,:,k),var(i,:,k),'o-','MarkerSize',4,'LineWidth',2,'Color',clut(i,:));
    end
%     for i =1:nStatus
%         try
%             fity(i,:) = interp1(1:nSteps1,peak(i,1:nSteps1,k),fitx,'pchip');
%             plot(fitx,fity(i,:),'-','LineWidth',3,'Color',clut(i,:));
%         end
%     end
    legend('D0_s','D0_r','D7_s','D7_r','Location','Best');
    legend('boxoff')
    axis tight
    ax3(k).XTick = xaxis;
    ax3(k).XTickLabel = names_dir;
    xlabel('Stimulus Orientation');
    ylabel('Normalized response amplitude')
    %     ylabel('Ratio(Day7/Day0)');
    title(['Preferred Orientation:' names_dir{k}  ',n=' num2str(count(k))])
end



ax3(nSteps+1) = subplot(2,nSteps1/2+1,nSteps1+2);hold on
%peak_dir size: 4*pref_dir
xaxis = [1:nSteps];% nOris1+1.5];
for i = 1:nStatus
    y = squeeze(peak(i,:,:)./peak(1,:,:));
    yr = squeeze(var(i,:,:)./peak(1,:,:));
    errorbar(xaxis+(i-2.5)*.1,peak_po(i,:),var_po(i,:),...
        'o','MarkerSize',8,'LineWidth',2,'Color',clut(i,:));
%     display(peak_dir(i,:)./min(peak_dir(1,:)));
end
% legend('D0_s','D0_r','D7_s','D7_r','Location','best');
% legend('boxoff')
axis tight
xlim([xaxis(1)-.5 xaxis(end)+.5])
ax3(nSteps+1).XTick = xaxis;
ax3(nSteps+1).XTickLabel = names_dir;
 xlim([xaxis(1)-.5 xaxis(end)+.5])
 
 linkaxes(ax3,'x')
 linkaxes(ax3,'y');
 
xlabel('Response@Preferred Orientation')
ylabel('Amplitude Ratio')
saveas(gcf,[figurename '_TC.fig'])
saveas(gcf,[figurename '_TC.png'])
