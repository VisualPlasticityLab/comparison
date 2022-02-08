function h0 = plotpeakselected_old(peak_selected)
ncell = size(peak_selected,3);
nSteps = size(peak_selected,2);
nSteps1 = nSteps - mod(nSteps,2);

[normal,pref_ori]=max(max(peak_selected));

for i=1:ncell
    resp1(:,:,i) = peak_selected(1:2,:,i)/normal(i);
    resp2(:,:,i) = peak_selected(3:4,:,i)/normal(i);
    
%      normal(i) = max(max(peak_selected(1:2,:,i)));
%     normal2(i) = max(max(peak_selected(3:4,:,i)));
%     resp1(:,:,i) = peak_selected(1:2,:,i)/normal(i);
%     resp2(:,:,i) = peak_selected(3:4,:,i)/normal2(i);
end

%%

thr = 0;
% resp (:,:,normal<thr) = NaN;
% resp (:,:,8) = NaN;
% resp (:,:,12) = NaN;

ori1 = (resp1(:,1:nSteps1/2,:)+resp1(:,nSteps1/2+1:nSteps1,:))/2;
ori2 = (resp2(:,1:nSteps1/2,:)+resp2(:,nSteps1/2+1:nSteps1,:))/2;
if nSteps> nSteps1
    ori1(:,nSteps1/2+1,:) = resp1(:,nSteps,:);
    ori2(:,nSteps1/2+1,:)= resp2(:,nSteps,:);
end
peak1 = nanmean(ori1,3);
peak2 = nanmean(ori2,3);

nOri = size(ori1,2);
%%
figurename=sprintf('%d cells with response greater than %.2f',sum(normal>thr),thr);
h0=figure('Position',[0 0 1200 1200],'Name',figurename);
ymin = min(resp1(:));
ymax = max(resp1(:));

h = (1:nOri)/nOri;
s = ones(1,nOri);
v = ones(1,nOri);
clut =hsv2rgb(h,s,v);

for i=1:nOri
if i<=nSteps1/2
    names{i} = sprintf('%.1fDegree',(i-1)/nSteps1*360);
else
    names{i} = 'blank';
end

subplot(nSteps1/4+1,2,i);hold on
for k=1:nSteps1/2
plot([k k+.4],peak1(:,i),'o-','MarkerSize',8,'LineWidth',3,'Color',clut(1,i,:))
plot([k+.1 k+.5],peak2(:,i),'x--','MarkerSize',8,'LineWidth',3,'Color',clut(1,i,:))
plot([k k+.4],squeeze(ori1(:,i,:)),'o-','Color', [.5 .5 .5 .5])
plot([k+.1 k+.5],squeeze(ori2(:,i,:)),'x-','Color', [.5 .5 .5 .5])
end

title(['Response:' names{i}])
legend('Day0','Day7','Location','best')
legend('boxoff')
axis tight
xlim([0.5 2.6])
end


subplot(nSteps1/4+1,2,nSteps1/2+2);hold on
for k=1:size(ori1,2)
plot([1 2],peak1(:,k),'o-','MarkerSize',8,'LineWidth',3,'Color',clut(1,k,:))
end
for k=1:size(ori1,2)
plot([1.1 2.1],peak2(:,k),'x--','MarkerSize',8,'LineWidth',3,'Color',clut(1,k,:))
end
legend(names,'Location','northwest')
legend('boxoff')
% axis tight
xlim([0.5 2.6])
