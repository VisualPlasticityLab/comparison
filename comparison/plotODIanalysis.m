function plotODIanalysis(fign,selected)

load([fign 'all.mat'],'PEAK0','PEAK00','PEAK2','ORI2','ODI0','ODI1',...
    'fig1img','F1cell','F1Contour','RGB1');

%% recalculate from eye1 and eye2 
if ~exist('PEAK0','var')
    load([fign 'all.mat'],'eye1','eye2');
    for jj=1:size(eye1.peak,3)-1 %the last one is background
        nth1= jj;
        nth2= jj;
        %F0 modulation
    %     eye1.peak(:,:,nth1)= eye1.peak(:,:,nth1)-nanmin(eye1.peak(:,1:end-1,nth1),[],2);
    %     eye1.peakS(:,:,nth1)=eye1.peakS(:,:,nth1)-nanmin(eye1.peakS(:,1:end-1,nth1),[],2);
    %     eye1.peakR(:,:,nth1)=eye1.peakR(:,:,nth1)-nanmin(eye1.peakR(:,1:end-1,nth1),[],2);
    %     
    %     eye2.peak(:,:,nth1)= eye2.peak(:,:,nth1)-nanmin(eye2.peak(:,1:end-1,nth1),[],2);
    %     eye2.peakS(:,:,nth1)=eye2.peakS(:,:,nth1)-nanmin(eye2.peakS(:,1:end-1,nth1),[],2);
    %     eye2.peakR(:,:,nth1)=eye2.peakR(:,:,nth1)-nanmin(eye2.peakR(:,1:end-1,nth1),[],2);    
    %     

    %     allpeaks = cat(1,eye1.peakR(:,:,nth1),eye2.peakR(:,:,nth2),eye1.peakS(:,:,nth1),eye2.peakS(:,:,nth2));
        [a,pref_ori]=nanmax(nanmax(cat(1,eye1.peak(:,:,nth1),eye2.peak(:,:,nth2))));
        [a1,ori1]=nanmax(eye1.peak(:,:,nth1));
        [a2,ori2]=nanmax(eye2.peak(:,:,nth2));

        [a_s,pref_ori_s]=nanmax(nanmax(cat(1,eye1.peakS(:,:,nth1),eye2.peakS(:,:,nth2))));
        [a_1s,pref_ori_1s]=nanmax(eye1.peakS(:,:,nth1));
        [a_2s,pref_ori_2s]=nanmax(eye2.peakS(:,:,nth2));

        [a_r,pref_ori_r]=nanmax(nanmax(cat(1,eye1.peakR(:,:,nth1),eye2.peakR(:,:,nth2))));
        [a_1r,pref_ori_1r]=nanmax(eye1.peakR(:,:,nth1));
        [a_2r,pref_ori_2r]=nanmax(eye2.peakR(:,:,nth2));

        [odi,~] = ODIcalculation(pref_ori,eye2.peak(:,:,nth2),eye1.peak(:,:,nth1));
        [odi_s,~] = ODIcalculation(pref_ori,eye2.peakS(:,:,nth2),eye1.peakS(:,:,nth1));
        [odi_r,~] = ODIcalculation(pref_ori,eye2.peakR(:,:,nth2),eye1.peakR(:,:,nth1));
        [odi_ss,~] = ODIcalculation(pref_ori_s,eye2.peakS(:,:,nth2),eye1.peakS(:,:,nth1));
        [odi_rr,~] = ODIcalculation(pref_ori_r,eye2.peakR(:,:,nth2),eye1.peakR(:,:,nth1));

        ORI0(jj) = pref_ori;
        ORI1(:,jj) = [ori1 ori2];
        ORI2(:,jj) = [pref_ori_s pref_ori_1s pref_ori_2s pref_ori_r pref_ori_1r pref_ori_2r];

        PEAK0(jj) =a;
        PEAK00(:,jj)= [eye1.peak(:,pref_ori,nth1) eye1.peakS(:,pref_ori,nth1) eye1.peakR(:,pref_ori,nth1)...
            eye2.peak(:,pref_ori,nth2) eye2.peakS(:,pref_ori,nth2) eye2.peakR(:,pref_ori,nth2)]; % use combined pref_ori during run/still 
        PEAK1(:,jj)= [a1 a2];
        PEAK2(:,jj)= [a_s a_1s a_2s a_r a_1r a_2r]; %use pref_ori_s and pref_ori_r 

        ST0(:,jj)= [eye1.error(:,pref_ori,nth1)  eye2.error(:,pref_ori,nth2)];
        ST1(:,jj)= [eye1.errorS(:,pref_ori,nth1) eye2.errorS(:,pref_ori,nth2)...
            eye1.errorR(:,pref_ori,nth1) eye2.errorR(:,pref_ori,nth2)];
        ST2(:,jj)= [eye1.errorS(:,pref_ori_s,nth1) eye2.errorS(:,pref_ori_s,nth2)...
            eye1.errorR(:,pref_ori_r,nth1) eye2.errorR(:,pref_ori_r,nth2)];

        ODI0(jj) = odi;
        ODI1(:,jj) = [odi_s  odi_r]; %using same pref_ori
        ODI2(:,jj) = [odi_ss  odi_rr]; %using run/still pref_ori

        try
            RGB1(:,jj) = ODI2rgb(odi,a); %ODIcolorbar to make a reference barmap
        end
    end
end
%% wanted to convert ODI to sqrt(ODI) but realized should do amp instead
% ODI0=convertODI(ODI0);
% ODI1=convertODI(ODI1);
%%
% selected = selected(PEAK0(selected)>.3);
segments = 7;
ncell = numel(selected);
numcell = num2str(numel(selected));

ori = round((max(ORI2(:))-1)/2);
colortype = ori;
h = (1:colortype)'/colortype;
s = ones(colortype,1);
v = ones(colortype,1);
clut =squeeze(hsv2rgb(h,s,v));
clut(ori+1:ori*2,:) = clut/2;
clut(ori*2+1,:) = [1 1 1]/2;
%clut2 =clut/2;

de=0:180/ori:360;
dx = cosd(de)*10;
dy = sind(de)*10;
dx(end) = 0; dy(end) = 0; % noise level sets to 0

level = 6;
peak = discretize(PEAK0,level)/level./PEAK0;

%% Plot ODI 
sets = 2; %define how many columns
h(2)=figure('Position',[200 200 300*sets 500]);
subplot(2,sets,1);hold on
plot([-1 1],[-1 1],'--','Color',[0.5 0.5 .5 .5])
scatter(ODI1(1,selected),ODI1(2,selected),abs(peak(selected)*50),'ko','filled'); 
scatter(nanmean(ODI1(1,selected)),nanmean(ODI1(2,selected)),100,'kx');
alpha(.3)

xlabel(['ODIstill \mu= ' sprintf('%.2f',nanmean(ODI1(1,selected)))])
ylabel(['ODIrun \mu= ' sprintf('%.2f',nanmean(ODI1(2,selected)))])
axis square
title('ODI@best direction:','FontSize',12)

subplot(2,sets,sets+1);hold on;
histplt2(ODI1(1,selected),ODI1(2,selected),-1.0:2/segments:1.0);
legend('still','run','Location','Northwest')
legend('Boxoff')
ylabel('Cell count')
title(['N=' numcell],'FontSize',12)
axis square

%PEAK2(:,jj)= [a_s a_1s a_2s a_r a_1r a_2r];
%%  plot seperate ODI
% subplot(2,sets,2);hold on;
% plot([-1 1],[-1 1],'--','Color',[0.5 0.5 .5 .5])
% scatter(ODI2(1,selected),ODI2(2,selected),peak(selected)*50,'ko','filled'); 
% scatter(nanmean(ODI2(1,selected)),nanmean(ODI2(2,selected)),100,'kx');
% alpha(.3)
% 
% xlabel(['ODIstill \mu= ',sprintf('%.2f',nanmean(ODI2(1,selected)))])
% ylabel(['ODIrun \mu= ',sprintf('%.2f',nanmean(ODI2(2,selected)))])
% axis square
% title('ODI@run/still direction','FontSize',12)
% 
% subplot(2,sets,2+sets);hold on;
% histplt2(ODI2(1,selected),ODI2(2,selected),-1.0:2/segments:1.0);
% legend('still','run','Location','Northwest')
% legend('Boxoff')
% ylabel('Cell count')
% title(['N=' numcell],'FontSize',12)
% axis square

%% plot raw data

maxpeak = max(prctile(PEAK2(:,selected),95));
subplot(2,sets,sets);hold on
% PEAK00(:,jj)= [eye1.peak(:,pref_ori,nth1) eye1.peakS(:,pref_ori,nth1) eye1.peakR(:,pref_ori,nth1)...
%         eye2.peak(:,pref_ori,nth2) eye2.peakS(:,pref_ori,nth2) eye2.peakR(:,pref_ori,nth2)]; 
hs = scatter(PEAK00(2,selected),PEAK00(5,selected),'bo');%still state contra-ipsi
hr = scatter(PEAK00(3,selected),PEAK00(6,selected),'ro');%run state contra-ipsi
hall = scatter(PEAK00(1,selected),PEAK00(4,selected),'ko','filled');%still state contra-ipsi
alpha(.3);
plot([0  maxpeak],[0 maxpeak],'--','Color',[0.5 0.5 .5 .5])
xlim([0  maxpeak])
xlabel('Ipsi eye response')
ylabel('Contra eye response')
axis square
title('Individual responses in ispi&contra-eyes','FontSize',12)

scatter(nanmean(PEAK00(2,selected),2),nanmean(PEAK00(3,selected),2),200,'bx');
scatter(nanmean(PEAK00(3,selected),2),nanmean(PEAK00(6,selected),2),200,'rx');
scatter(nanmean(PEAK00(1,selected),2),nanmean(PEAK00(4,selected),2),200,'kx');
legend({'still','run','overall'})
legend('boxoff')

subplot(2,sets,sets*2);hold on
histogram(ODI0(selected),-1.0:2/segments:1.0);
legend('overall','Location','Northwest')
legend('Boxoff')
ylabel('Cell count')
title(['ODI \mu= ',sprintf('%.2f',nanmean(ODI0(selected)))])
% maxpeak = max(prctile(PEAK2(:,selected),95));
% subplot(2,sets,3);hold on
% plot([0  maxpeak],[0 maxpeak],'b--')
% plot([0 -maxpeak],[0 maxpeak],'b--')
% %PEAK2(:,jj)= [a_s a_1s a_2s a_r a_1r a_2r];
% for jj=1:numel(selected) % plot orientation for still and run
%     plot([-PEAK2(2,jj) PEAK2(3,jj)],[PEAK2(5,jj) PEAK2(6,jj)],'.-','Color',[.5 .5 .5]);
% end
% plot([-mean(PEAK2(2,selected),2) mean(PEAK2(3,selected),2)],...
%     [mean(PEAK2(5,selected),2) mean(PEAK2(6,jj))],'o-','Color','k');
% 
% ylabel('Response run')
% ylim([0 maxpeak])
% %set(gca,'YTickLabel','')
% set(gca,'YAxisLocation','origin'); 
% xlabel('Response still')
% set(gca,'XTick',[-.25  .25])
% set(gca,'XTickLabel',{'Ispi','Contra'})
% title('Averaged responses for contra/ispi-eye during run/still')
% 
% subplot(2,sets,sets+3);
% histplt2([-PEAK2(2,selected) PEAK2(3,selected)],[ -PEAK2(5,selected) PEAK2(6,selected)],[-.05:1/200:.05]);
% set(gca,'YAxisLocation','origin'); 
% legend('still','run')
% legend('boxoff')
% ylabel('Response run')
% xlabel('Response still')
% title('Preferred orientation:run/still','FontSize',16)


% subplot(2,3,3);hold on;
% plot([-1 1],[-1 1],'--')
% scatter(ODI2(1,selected),ODI2(2,selected),peak(selected)*20,'*'); 
% scatter(nanmean(ODI2(1,selected)),nanmean(ODI2(2,selected)),200,'*');
% ylabel('ODIrun')
% xlabel('ODIstill')
% axis square
% title('Preferred orientation:run/still')
% 
% subplot(2,3,6);hold on;
% histplt2(ODI2(1,selected),ODI2(2,selected),-1.0:2/7:1.0);
% legend('ODIstill','ODIrun')
% legend('Boxoff')
% ylabel('Cell count')
% title(['N=' numcell],'FontSize',16)
% axis square
%%
h(1)=figure('Position',[200 600 800 600]);

subplot('Position',[0.15 0.05 0.8 .05]);hold on
odi_bar= repmat(-1:.005:1,21,1);
peak_column = repmat((0:.05:1)',1,401);
ODIbar = ODI2rgb(odi_bar,peak_column);
imshow(ODIbar,'InitialMagnification',400)

subplot('Position',[0.1 0.2 0.8 .7])
SIZ=size(fig1img);hold on;
x=[];
y=[];
sz0=[];
sz1=[];
sz2=[];
clut0=[];
clut1=[];
clut2=[];
dir1=[];
dir2=[];

for jj=selected
    x(end+1)= F1cell(jj).Position(1);
    y(end+1)= F1cell(jj).Position(2);
    sz0(end+1) = PEAK0(jj).*peak(jj)*15;%peak at all
    clut0(end+1,:) = ODI2rgb(ODI0(jj),sz0(end));% all ODI
    sz1(end+1) = PEAK2(1,jj).*peak(jj)*15;%peak at still
    dir1(end+1) = ORI2(1,jj);
    clut1(end+1,:) = ODI2rgb(ODI1(1,jj),sz1(end));%still ODI
    sz2(end+1) = PEAK2(4,jj).*peak(jj)*15;%peak at run
    dir2(end+1) = ORI2(4,jj);
    clut2(end+1,:) = ODI2rgb(ODI1(2,jj),sz2(end)); %run ODI
end
hh1=scatter(x,y,100*sz0,clut0,'filled'); 
% alpha(hh1,.5);
sz1(find(sz1==0)) = 0.001;
sz2(find(sz2==0)) = 0.001;
scatter(x,y,100*sz1,clut1,'LineWidth',2)
scatter(x,y,100*sz2,clut2,'LineWidth',3);
axis off;axis tight;axis ij

hh = gca; % plot orientation legend
x0 = hh.XLim(1);
y0 = hh.YLim(1);
oricolormap(x0+10*ori,y0+10*ori,16*ori,ori,'still','run')
alpha(.3)

% for jj=1:numel(selected) % plot orientation for still and run
%     arrow([x(jj) y(jj)],[x(jj)+dx(dir1(jj))*sqrt(sz1(jj)) y(jj)+dy(dir1(jj))*sqrt(sz1(jj))],20,...
%         'EdgeColor',[0 0 0],'FaceColor',[1 1 1],'LineWidth',1,'Length',8,'tipangle',30);
%     %'EdgeColor',clut(dir1(jj),:)
%     arrow([x(jj) y(jj)],[x(jj)+dx(dir2(jj))*sqrt(sz2(jj)) y(jj)+dy(dir2(jj))*sqrt(sz2(jj))],20,...
%         'EdgeColor',[0 0 0],'FaceColor',[.5 .5 .5],'LineWidth',1,'Length',8,'tipangle',30);
% %'EdgeColor',[0 0 0],'FaceColor',clut(dir2(jj),:)
% end
axis off
scalebar('Unit','pixel')
title(['N=' numcell],'FontSize',16)
% close(h(1))
%% save files
pngname = [fign 'selected' numcell];
if exist([pngname '.png'],'file')
    pngname = [pngname '_'];
end
        
saveas(h(1),[pngname '.png']);
saveas(h(2),[pngname '_ODI.png']);

h(3)=figure;
imshow(fig1img,[]);hold on;
for jj=selected
    contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,jj),SIZ(1),SIZ(2)),[0.01 1],'LineColor',RGB1(:,jj),'LineWidth',1)
%     text(F1cell(jj).Position(1),F1cell(jj).Position(2),F1cell(jj).String,...
%         'Color',[1 0 0],'FontSize',16)
end
imshow(ODIbar,'InitialMagnification',400)
saveas(h(3),[pngname '_cells.png']);

