%%%%calculate_OSI: adapted from calculate_ODI5: not considering run/still but first plot
%%%%the averaged traces for all conditions 10-08-2018

clear all;close all;clc
%% load figures and select a ROI to print and analyze
[figname1,path1a]=uigetfile('.fig','find the fig and the memmap.mat');
fign0=strtok(figname1,'t');

figimg = sbxshownplane;
if menu('Make sure 1st recording file is showed?','Yes','reselect') ==2
    figimg = sbxshownplane;
end
if iscell(figimg)
    fig1img = figimg{menu('which plane to select?','1st','2nd','3rd')};
else
    fig1img = figimg;
end
close(gcf)
%%
try
    fig1=openfig(fullfile(path1a,figname1));
catch
    [day1sig,p1]=uigetfile('.signals','select signal file for day1eye');
    s1 = matfile(fullfile(p1,day1sig));
    %s1.A has all the contours and text position is   com(Aor(:,1:end),d1,d2);
    c1 = plot_contours(s1.A,s1.V,0.95,1); %thr=.95,display=TRUE
    fig1=gcf;
%     saveas(fig1,[figname1 '.fig'],'fig');
%     saveas(fig1,[figname1 '.png'],'png');
end

%%
c = calcontour(F1Contour);
[f,p]=uigetfile('*.mat','Load region from previous file?');
if f~= 0
    load(fullfile(p,f),'maskpic0');
else
    figure(fig1);
    maskpic0=roipoly;
end
close(fig1)

[prefilename,path]=uigetfile('*.mat','Pick the pre-matched pairs?');
if f~= 0
    load(fullfile(p,f),'maskpic0','pair1','pair2','F1cell','F2cell',...
    'F1Contour','F2Contour');
end

% Set potential good cells criteria: response and size
try 
    memmapfile=matfile(fullfile(path1a,[fign0(1:pos(end-1)) 'memmap.mat']));
    magnification = memmapfile.magnification;
catch
    magnification =2;
end
if mod(magnification,1)|| magnification==2
    magnification = magnification*2.5;
end 
load('magnificationlist');
mag= maglist(magnification)

c.area_cutoff = mag^2*30 *[.55 3.3] ; % cellsize = [96 396] under 2X
c.circ_cutoff = 1.75; % 1/3 < a/b < 3 c.circ=5/3
c.skew_cutoff = 1; 
  
pair0 = reshape(maskpic0,1,prod(SIZ))*F1Contour >0;
pair1 = c.area<=c.area_cutoff(2) & c.area>=c.area_cutoff(1) & c.circ<=c.circ_cutoff;
% responselimit = .03;
% pair1 = (PEAK0>=responselimit)&(SNR(1:end-1)>1.3)&pair1;
%%
h21=figure;
imshow(fig1img)%ipsi eye,red;contra eye,green
hold on;
for jj=1:ncell
    contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,jj),SIZ(1),SIZ(2)),[0.05 1],'LineColor',[pair0(jj) pair0(jj)&~pair1(jj) 0 ],'LineWidth',2);
%      text(F1cell(jj).Position(1),F1cell(jj).Position(2),[F1cell(jj).String sprintf('=%d,%.1f',c.area(jj),c.circ(jj))],...
%         'Color',[pair0(jj)&pair1(jj) 0 0],'FontSize',16)
end
saveas(h21,'cellshape.fig')
saveas(h21,'cellshape.png')

%% load sig data
[eye1,eye2,sigtype]=loadsigdata; %sigtype: ca, sp, ML
plotrunrep(eye1.matrix,eye2.matrix)
SNR = cal_SNR;%%

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
    
%%
fign = [fign0 sigtype '_'];
save([fign 'all.mat'],'SNR','eye1','eye2','ORI0','ORI1','ORI2','PEAK0','PEAK00','PEAK1','PEAK2','ST0','ST1','ST2','ODI0','ODI1','ODI2',...
    'fig1img','F1cell','F1Contour','SIZ','RGB1','c','pair1','pair0','maskpic0','fign','-v7.3');
%%
selected1 = find(pair0&pair1);
% Analysis Summary for automatic selection
plotODIanalysis(fign,selected1)
% close all
 plotTC(fign,selected1);
%ManuelCellSelection
