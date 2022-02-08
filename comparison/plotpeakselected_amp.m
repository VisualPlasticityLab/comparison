function h0 = plotpeakselected_amp(peak_selected)
%peqk_selected: status(4:day1still,day1run,day2still,day2run)*nSteps*ncell

ncell = size(peak_selected,3);
nSteps = size(peak_selected,2);
nSteps1 = nSteps - mod(nSteps,2);
nOris1 = nSteps1/2;
nOris = nOris1+(nSteps-nSteps1);
% nOris1 = nSteps1;edi
% nOris = nSteps;

% [max_amp,pref_ori0]=nanmax(peak_selected(3,:,:));
[max_amp,pref_ori0]=nanmax(nanmax(peak_selected(:,:,:)));
pref_ori = pref_ori0;
pref_ori(pref_ori>nOris1) = pref_ori(pref_ori>nOris1)-nOris1;

ori = squeeze(nanmean(reshape(peak_selected(:,1:nSteps1,:),[],nOris1,2,ncell),3));
% ori = peak_selected;
if nOris
    ori(:,nOris,:) = peak_selected(:,nSteps,:);
end
% [~,pref_ori] = max(max(ori));

for i=1:ncell
    amp(i) = ori(1,pref_ori(i),i);
    resp(:,:,i) = ori(:,:,i)./amp(i);
%     amp(i) = ori(1,pref_ori(i),i);
%     resp(1,:,i) = real(log2(ori(3,:,i)./amp(i)));
%     amp2(i) = ori(3,pref_ori(i),i);
%     resp(2,:,i) = real(log2(ori(4,:,i)./amp2(i)));
    
end
%% Only max_amp > thr and not isinf/isnan in any days
thr = 0.05;
% resp (:,:,normal<thr) = NaN;
resp(:,:,max_amp<thr) = NaN;
resp(isinf(resp)) = NaN;
% thr_up = prctile( resp(:),[1 99]);
% resp(resp>thr_up(2)) =NaN;
% resp(resp<thr_up(1)) =NaN;
% resp(:,:,isnan(sum(sum(resp,1)))) = NaN;

for i=1:nOris
    selected=(pref_ori==i)&(~isnan(sum(sum(resp))));
    count(i) = sum(selected);
    peak(:,:,i) = nanmean(resp(:,:,selected),3);
    var(:,:,i) = nanstd(resp(:,:,selected),0,3)/sqrt(count(i));
end

%%
h = (1:nOris)/nOris;
s = ones(1,nOris);
v = ones(1,nOris);
clut =hsv2rgb(h,s,v)/2;

for i=1:nOris1
    if mod((i-1)*360,nSteps1)~=0
        names{i} = sprintf('%.1f/%.1f',(i-1)*360/nSteps1,(i-1)*360/nSteps1+180);
    else
        names{i} = sprintf('%d/%d',(i-1)*360/nSteps1,(i-1)*360/nSteps1+180);
    end
end
if nOris>nOris1
    names{nOris} = 'blk';
end

%% plot Response Amplitude
ymax = max(peak(:));

figurename=sprintf('Response Amplitude:%d cells with response greater than %.2f',sum(count),thr);
h0=figure('Position',[200 200 1200 600],'Name',figurename);
xaxis=linspace(0,.5,size(peak,1)); %[k;k+.4];
for i=1:nOris
    subplot(nOris1/2+1,2,i);hold on
    for k=1:nOris %pref_ori
        errorbar(xaxis+k,squeeze(peak(:,i,k)),squeeze(var(:,i,k)),'o-','MarkerSize',4,'LineWidth',2,'Color',clut(1,i,:))
        try
            plot(xaxis+k,squeeze(resp(:,i,pref_ori==k)),'o-','LineWidth',.3,'Color', [.5 .5 .5 .5])
        end
        errorbar(xaxis+k,squeeze(peak(:,i,k)),squeeze(var(:,i,k)),'o-','MarkerSize',4,'LineWidth',2,'Color',clut(1,i,:))
        [T(k,i),P(k,i)]=ttest(resp(1,i,pref_ori==k));
    end
    
    legend('Day0--Day7','Location','northeast')
    legend('boxoff')
%     axis tight
    
    ax=gca;
    % when use log2
    ax.YScale = 'linear'; 
%     ylim([-5 5]);
%     ax.YTick = -4:2:4;
%     ax.YTickLabel = 2.^ax.YTick;
    %
    % ax.YScale = 'log';
    % ylim([.025 40]);
    % ax.YTick = 10.^(-1.6:1.6:1.6);
    % ax.YTickLabel = [.025,1, 40];
    
    xlim([0.5 nOris+1])
    ylim( [-1 ceil(ymax)])

    ax.XTick = (1:nOris)+.25;
    ax.XTickLabel = names;
    
    xlabel('Preferred Orientation');
    ylabel('Ratio(Day7/Day0)');
    title(['Response@' names{i}])
    
end

disp(T)
disp(P)

subplot(nOris1/2+1,2,nOris1+2);hold on
% figure;hold on
for k=1:nOris
    errorbar(xaxis+k,squeeze(peak(:,k,k)),squeeze(var(:,k,k)),'o-','MarkerSize',4,'LineWidth',2,'Color',clut(1,k,:))
end
xlim([0.5 nOris+1])
ylim( [-1 ceil(ymax)])
xlabel('Preferred Orientation')
ylabel('Amplitude Ratio')
ax=gca;
ax.XTick = 1:nOris;
ax.XTickLabel = names;

