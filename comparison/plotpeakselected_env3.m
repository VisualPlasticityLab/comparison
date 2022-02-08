function [h0,ori_ba,resp] = plotpeakselected_env3(peak_selected)
%% peak_selected: status*nSteps*ncell
% status(2or4:Post-still,Post-run,day2still,day2run)

%% discard cells with the max response <threshold
thr = 0.00;
% peak_selected(peak_selected<0) = 0;
peak_selected(:,:,max(max(peak_selected))<thr) = [];

%%
nStatus = size(peak_selected,1);
nSteps = size(peak_selected,2);
nSteps1 = nSteps - mod(nSteps,2);
nOris1 = nSteps1/2;
nOris = nOris1+(nSteps-nSteps1);

[amp_max1,~]=nanmax(nanmax(peak_selected(3:4,:,:)));
[amp_min1,~]=nanmin(nanmin(peak_selected(3:4,:,:)));
[~,pref_ori1]=nanmax(nanmax(peak_selected(3:4,:,:)));
pref_ori1(pref_ori1>nOris1) = pref_ori1(pref_ori1>nOris1)-nOris1;

[amp_max0,~]=nanmax(nanmax(peak_selected(1:2,:,:)));
[amp_min0,~]=nanmin(nanmin(peak_selected(1:2,:,:)));
[~,pref_ori0]=nanmax(nanmax(peak_selected(1:2,:,:)));
pref_ori0(pref_ori0>nOris1) = pref_ori0(pref_ori0>nOris1)-nOris1;

%% use Pre- or Post- as the preferred orientation
% pref_ori = pref_ori1;
% amp_max = amp_max1(:);
% amp_min = amp_min1(:);
[amp_max,~]=nanmax(nanmax(peak_selected));
[amp_min,~]=nanmin(nanmin(peak_selected));
%% Only not isinf/isnan in any days
%resp size: 2/4*nSteps*cells
%peak size: 2/4*nSteps*pref_ori
%peak_ori size: 2/4*pref_ori

for i=1:size(peak_selected,3)
    %     resp(:,:,i) = (peak_selected(:,:,i)-amp_min(i))/(amp_max(i)-amp_min(i));
    %     resp(:,:,i) = peak_selected(:,:,i);
    resp(1:2,:,i) = (peak_selected(1:2,:,i)-amp_min0(i))/(amp_max0(i)-amp_min0(i));
    resp(3:4,:,i) = (peak_selected(3:4,:,i)-amp_min1(i))/(amp_max1(i)-amp_min1(i));
end
% resp = real(log2(resp+1));
resp(isinf(resp)) = NaN;
% resp_ori = squeeze(nanmean(reshape(resp(:,1:nSteps1,:),4,nOris1,2,ncell),3));
% if nOris>nOris1
%     resp_ori(:,nOris,:) = resp(:,nSteps,:);
% end

%%
colortype = 2; % before &after
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
%% plot Overall changes

for b=1:nOris
    for a=1:nOris
        ori_ba(b,a)=sum(pref_ori0==b&pref_ori1==a);
    end
end
oria = repmat(1:nOris,nOris,1);
orib = oria';
for k = 1:nOris
    ori_b(k) =sum(pref_ori0==k);
    ori_a(k) =sum(pref_ori1==k);
end
%%
% imagesc(ori_ba/repmat(ori_b,[
xaxis= 1:nSteps;% nSteps1+2.5];

figure('Position',[390 100 370 322]);hold on
subplot(1,2,1)
runstim = 45/(180/nOris1)+1;
stillstim = 135/(180/nOris1)+1;
nonexp = find(oria~=runstim&oria~=stillstim&ori_ba~=0);
runexp = find(oria==runstim&ori_ba~=0);
stiexp = find(oria==stillstim&ori_ba~=0);
scatter(orib(nonexp),oria(nonexp),ori_ba(nonexp)*15,'filled','MarkerFaceColor',[.8 .8 .8]);
scatter(orib(runexp),oria(runexp),ori_ba(runexp)*15,'filled','MarkerFaceColor',[.5 0 0]+.5);
scatter(orib(stiexp),oria(stiexp),ori_ba(stiexp)*15,'filled','MarkerFaceColor',[0 0 .5]+.5);
for b=1:nOris
    for a=1:nOris
        text(b,a,sprintf('%d%%',round(ori_ba(b,a)/ori_b(b)*100)),'HorizontalAlignment','Center');
        %     text(b,a,sprintf('%d%',ori_ba(b,a)),'HorizontalAlignment','Center');
    end
end
axis([0.5 nOris+.5 0.5 nOris+.5])
bax =gca;
bax.XTick = xaxis;
bax.XTickLabel = names_ori;
xlabel('Pre-exposure best ori(degree)')
bax.YTick = xaxis;
bax.YTickLabel = names_ori;
ylabel('Post-exposure best ori(degree)')
subplot(1,2,2)
for i = 1:nStatus
    errorbar(xaxis,nanmean(resp(i,:,:),3),nanstd(resp(i,:,:),0,3)/sqrt(size(resp,3))...
        ,'o-','MarkerSize',4,'LineWidth',2,'Color',clut(1,i,:))
    
    %         try
    %             fity(i,:) = interp1(1:nSteps1,nanmean(resp(i,1:nSteps1,cells),3),fitx,'pchip');
    %             plot(fitx,fity(i,:),'-','LineWidth',1,'Color',clut(1,i,:));
    %         end
end
xlim([xaxis(1)-.5 xaxis(end)+.5])


%% plot Tuning curve
%resp size: 4*nSteps*cells
%peak size: 4*nSteps*pref_ori
figurename=sprintf('Tuning curve_%d cells with response greater than %.2f',size(peak_selected,3),thr);
h0=figure('Position',[100 100 800 600],'Name',figurename);
fitx = linspace(1,nSteps1,nSteps1*3);
[xpos,ypos,xwidth,yheight] = figurepara(nOris,nOris);
for k=1:nOris^2
    hn(k)=subplot('Position',[xpos(mod(k-1,nOris)+1),ypos(ceil(k/nOris)),xwidth,yheight]); % plot tuning curve @ preferred orientation
    hold on
    cells = (pref_ori0==orib(k)&pref_ori1==oria(k));
    for i = 1:nStatus
        errorbar(xaxis,nanmean(resp(i,:,cells),3),nanstd(resp(i,:,cells),0,3)/sqrt(sum(cells))...
            ,'o-','MarkerSize',4,'LineWidth',2,'Color',clut(1,i,:))
        %         try
        %             fity(i,:) = interp1(1:nSteps1,nanmean(resp(i,1:nSteps1,cells),3),fitx,'pchip');
        %             plot(fitx,fity(i,:),'-','LineWidth',1,'Color',clut(1,i,:));
        %         end
    end
    text(nOris,-.05,sprintf('n=%d(%d%%)',sum(cells)),round(ori_ba(b,a)/ori_b(b)*100,'HorizontalAlignment','Center')
    axis tight
    ax=gca;
    xlim([xaxis(1)-.5 xaxis(end)+.5])
    axis off
    %     ax.XTick = xaxis;
    %     ax.XTickLabel = names;
    %     if k==nOris1/2+1
    %         legend({'pre-S','pre-R','post-S','post-R'},'Orientation','Horizontal','Location','best');
    %         legend('boxoff')
    %     end
    %     xlabel('Stimulus Orientation');
    %     ylabel('Normalized response')
    %     title(['Preferred Orientation:' names{k} '&' names{k+nOris1} 'n=' num2str(count(k))])
end
linkaxes(hn,'y')
%
% saveas(gcf,[figurename '.fig'])
% saveas(gcf,[figurename '.png'])
