function h0 = plotpeakselected_env(peak_selected)
%% peak_selected: status(2or4:Post-still,Post-run,day2still,day2run)*nSteps*ncell
thr = 0.05;

peak_selected(peak_selected<0) = 0;
peak_selected(:,:,max(max(peak_selected))<thr) = [];


nStatus = size(peak_selected,1);
nSteps = size(peak_selected,2);
nSteps1 = nSteps - mod(nSteps,2);
nOris1 = nSteps1/2;
nOris = nOris1+(nSteps-nSteps1);
% nOris1 = nSteps1;
% nOris = nSteps;
%ori = squeeze(nanmean(reshape(peak_selected(:,1:nSteps1,:),4,nOris1,2,ncell),3));
% if nOris>nOris1
%     ori(:,nOris,:) = peak_selected(:,nSteps,:);
% end
% [~,pref_ori] = max(max(ori));
% for i=1:ncell
%     amp_max(i) = max_amp(i);
%     amp_min(i) = min(ori(1,:,i));
%     resp(:,:,i) = (ori(:,:,i)-amp_min(i))./(amp(i)-amp_min(i));
% end
% resp = real(log2(resp));
[amp_max1,pref_ori1]=nanmax(nanmax(peak_selected(3:4,:,:)));
[amp_min1,~]=nanmin(nanmin(peak_selected(3:4,:,:)));
[amp_max0,pref_ori0]=nanmax(nanmax(peak_selected(1:2,:,:)));
[amp_min0,~]=nanmin(nanmin(peak_selected(1:2,:,:)));
%  [amp_max0r,~]=nanmax(peak_selected(2,:,:));
%  [amp_min0r,~]=nanmin(peak_selected(2,:,:));
pref_ori1(pref_ori1>nOris1) = pref_ori1(pref_ori1>nOris1)-nOris1;
pref_ori0(pref_ori0>nOris1) = pref_ori0(pref_ori0>nOris1)-nOris1;

 % use Pre- or Post- as the preferred orientation
pref_ori = pref_ori0;
amp_max = amp_max0(:);
amp_min = amp_min0(:);

% th_max=prctile(amp_max(:),[2 98]);
% resp(:,:,amp_max<th_max(1)|amp_max>th_max(2)) = [];
% 
% th_min=prctile(amp_min,[2 98]);
% amp_min(amp_min<th_min(1)|amp_min>th_min(2))= NaN;

%% Only max_amp > thr and not isinf/isnan in any days
for i=1:size(peak_selected,3)
    resp(:,:,i) = (peak_selected(:,:,i)-amp_min(i))/(amp_max(i)-amp_min(i));
%    resp(:,:,i) = peak_selected(:,:,i);
end
% resp = real(log2(resp+1));

% resp (:,:,normal<thr) = NaN;
resp(isinf(resp)) = NaN;
% % thr_up = prctile(resp(:),[1 99]);
% % resp(resp<thr_up(1)|resp>thr_up(2)) =NaN;
% resp(:,:,isnan(sum(sum(resp,1)))) = NaN;

% resp_ori = squeeze(nanmean(reshape(resp(:,1:nSteps1,:),4,nOris1,2,ncell),3));
% if nOris>nOris1
%     resp_ori(:,nOris,:) = resp(:,nSteps,:);
% end

%resp size: 4*nSteps*cells
%peak size: 4*nSteps*pref_ori
%peak_ori size: 4*pref_ori

%%
% colortype = nOris;
colortype = 2;
h = (1:colortype)/colortype;
s = ones(1,colortype);
v = ones(1,colortype);
clut =hsv2rgb(h,s,v);
clut(:,3:4,:) =hsv2rgb(h,s,v)/2;

for i=1:nSteps1
    if mod((i-1)*360,nSteps1)==0
        names{i} = sprintf('%d',(i-1)*360/nSteps1);
    else
        names{i} = sprintf('%.1f',(i-1)*360/nSteps1);
    end
end
if nOris>nOris1
    names{end+1} = 'blk';
end

for i=1:nOris1
    if mod((i-1)*360,nSteps1)~=0
        names_ori{i} = sprintf('%.1f/%.1f',(i-1)*360/nSteps1,(i-1)*360/nSteps1+180);
    else
        names_ori{i} = sprintf('%d/%d',(i-1)*360/nSteps1,(i-1)*360/nSteps1+180);
    end
end
if nOris>nOris1
    names_ori{nOris} = 'blk';
end

for i=1:nOris
    selected{i}=squeeze((pref_ori==i));%&(~isnan(sum(sum(resp))));
    count(i) = sum(selected{i});
    if i<nOris1+1
        pref_resp{i} = reshape(resp(:,[i i+nOris1],selected{i}),[],2*count(i));
    else
        pref_resp{i} = reshape(resp(:,i+nOris1,selected{i}),[],count(i));
    end
    peak_ori(:,i)= nanmean(pref_resp{i},2);
    var_ori(:,i) = nanstd(pref_resp{i},1,2)/sqrt(count(i));
    peak(:,:,i) = nanmean(resp(:,:,selected{i}),3);
    var(:,:,i) = nanstd(resp(:,:,selected{i}),0,3)/sqrt(count(i));
end
%% plot Overall changes
figure('Name','Population Histogram','Position',[600 400 1000 550]);
subplot(2,3,1);%plot(squeeze(pref_ori0),squeeze(pref_ori1),'o');
histplt2(squeeze(pref_ori0),squeeze(pref_ori1));
ax=gca;
ax.XTick = 1:nSteps;
ax.XTickLabel = names_ori;
ylabel('Cell count')
title('preferred orientation');
legend('Pre-','Post-');
legend('boxoff')
% subplot(2,2,2);histogram(pref_ori1);title('preferred direction')
subplot(2,3,2);histplt2(amp_max0,amp_max1,[0:max(amp_max1)/10:max(amp_max1)]);
title('amp\_max');legend('Pre-','Post-');legend('boxoff')
subplot(2,3,3);histplt2(amp_min0,amp_min1,[0:max(amp_min1)/10:max(amp_min1)]);
title('amp\_min');legend('Pre-','Post-');legend('boxoff')
for k=1:nOris
    ax(k) = subplot(2,nOris1+1,k+nOris1+1);hold on
%     histplt2(resp(1,k,selected{k}),resp(2,k,selected{k}),linspace(0,1,10));
    plot([1;2], squeeze([resp(1,k,selected{k}) ; resp(2,k,selected{k})] ),'x-','Color',[.5 .5 .5 .3],'LineWidth',.5);
    plot([1;2], [nanmean(resp(1,k,selected{k})) ; nanmean(resp(2,k,selected{k}))],'ko-','LineWidth',1);
    xlim([.5 2.5])
    set(gca,'XTick',[1 2])
    set(gca,'XTickLabel',{'Pre-','Post-'})
    xlabel(names_ori(k))
    ylabel('Maxium response (df/F)')
end
linkaxes(ax,'y')
ylim([-.05 1.5])
saveas(gcf,'Overall_stats.fig')

%% plot Tuning curve
%resp size: 4*nSteps*cells
%peak size: 4*nSteps*pref_ori
figurename=sprintf('Tuning curve_%d cells with response greater than %.2f',sum(count),thr);
h0=figure('Position',[100 100 1000 550],'Name',figurename);
xaxis = [1:nSteps];% nSteps1+2.5];
fitx = linspace(1,nSteps1,nSteps1*3);


for k=1:nOris
    hn(k)=subplot(2,nOris1/2+1,k); % plot tuning curve @ preferred orientation
    hold on
    for i = 1:nStatus
        errorbar(xaxis,peak(i,:,k),var(i,:,k),'o-','MarkerSize',8,'LineWidth',2,'Color',clut(1,i,:))
    end

    for i =1:nStatus
        %         fity(i,:) =  fit2peakgaussion(k,peak0)
        %         fity(i,:) = interp1(xaxis,peak(i,1:nSteps1,k),fitx,'linear');
        %         plot(xaxis,squeeze(resp(i,:,selected{k})),'x','MarkerSize',3,'LineWidth',.1,'Color', [squeeze(clut(1,i,:));.25])
        try
            fity(i,:) = interp1(1:nSteps1,peak(i,1:nSteps1,k),fitx,'pchip');
            %             plot(fitx,fity(i,:),'-','LineWidth',3,'Color',clut(1,i,:));
        end
    end
    %     [T(k,i),P(k,i)]=ttest(resp(1,i,pref_ori==i));
    %     axis tight
    ax=gca;
    %     % when use log2
    %     ax.YScale = 'linear';
    %     ylim([-8 4]);
    %     ax.YTick = -8:2:4;
    %     ax.YTickLabel = 2.^ax.YTick;
    %     % when use real ratio
    %         ax.YScale = 'log';
    %         ylim([1/16 16]);
    %         ax.YTick = 2.^(-4:2:4);
    %         ax.YTickLabel = 2.^(-4:2:4);
    xlim([xaxis(1)-.5 xaxis(end)+.5])
%     ylim([-.25 1.25])
    ax.XTick = xaxis;
    ax.XTickLabel = names;
    if k==nOris1/2+1
            legend({'pre-S','pre-R','post-S','post-R'},'Orientation','Horizontal','Location','best');
        legend('boxoff')
    end
    xlabel('Stimulus Orientation');
    ylabel('Normalized response')
    %     ylabel('Ratio(Day7/Pre-)');
    title(['Preferred Orientation:' names{k} '&' names{k+nOris1} 'n=' num2str(count(k))])
end

% disp(T)
% disp(P)
linkaxes(hn,'y')


subplot(2,nOris1/2+1,nOris1+2);hold on
% for k=1:nSteps
%     m=k-(k>nOris1)*nOris1;
%     for i=1:4
%     errorbar(xaxis(k)+(i-2.5)*.1,squeeze(peak(i,k,m)),squeeze(var(i,k,m)),'o','MarkerSize',4,'LineWidth',2,'Color',clut(1,i,:))
%     end
% end

%resp_ori size: 4*cells
%peak_ori size: 4*pref_ori

xaxis = [1:nOris];% nOris1+1.5];
for i = 1:nStatus
    y = squeeze(peak(i,:,:)./peak(1,:,:));
    yr = squeeze(var(i,:,:)./peak(1,:,:));
%     errorbar(xaxis+(i-2.5)*.1,y,yr,'o','MarkerSize',8,'LineWidth',2,'Color',clut(1,i,:));

%     errorbar(xaxis+(i-2.5)*.1,peak_ori(i,:)./min(peak_ori(1:2,:)),var_ori(i,:),'o','MarkerSize',8,'LineWidth',2,'Color',clut(1,i,:));
    errorbar(xaxis+(i-2.5)*.1,peak_ori(i,:),var_ori(i,:),'o','MarkerSize',8,'LineWidth',2,'Color',clut(1,i,:));
end
% legend('D0_s','D0_r','D7_s','D7_r','Location','best');
% legend('boxoff')
ax=gca;
% %when use log2
% ax.YScale = 'linear';
% ylim([-4 2]);
% ax.YTick = -4:2:2;
% ax.YTickLabel = 2.^ax.YTick;

%when use real ratio
%     ax.YScale = 'log';
%     ylim([1/8 8]);
%     ax.YTick = 2.^(-3:2:3);
%     ax.YTickLabel = 2.^(-3:2:3);
xlim([xaxis(1)-.5 xaxis(end)+.5])
% ylim([-.25 1.25])
ax.XTick = xaxis;
ax.XTickLabel = names_ori;
xlabel('Response@Preferred Orientation')
ylabel('Amplitude Ratio')
% saveas(gcf,[figurename '.fig'])
% saveas(gcf,[figurename '.png'])
