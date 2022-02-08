%%%%calculate_ODI5
%%%%update from calculate_ODI4: not considering run/still but first plot
%%%%the averaged traces for all conditions 08-07-2018
clear all;close all;clc
%% load figures and select a ROI to print and analyze
[figname1,path1a]=uigetfile('.fig','find the fig and the memmap.mat');
if ~figname1
    [figname1,path1a]=uigetfile('.fig','find the *.signals instead');
end
fign0=strtok(strtok(figname1,'.'),'t');

% try
%     figimg = sbxshownplane;
%     if menu('Make sure 1st recording file is showed?','Yes','reselect') ==2
%         figimg = sbxshownplane;
%     end
%     if iscell(figimg)
%         fig1img = figimg{menu('which plane to select?','1st','2nd','3rd')};
%     else
%         fig1img = figimg;
%     end
% catch
    h=openfig(fullfile(path1a,figname1))
    fig1img =h.Children.Children(end).CData;
    fig1img = fig1img/prctile(fig1img(:),97);
% end
%close(gcf)
%% gather each cell's contour
fig1=h;
% try
%     fig1=openfig(fullfile(path1a,figname1));
% catch
%     [day1sig,p1]=uigetfile('.signals','cannot open figure, select *.signals instead');
%     s1 = matfile(fullfile(p1,day1sig));
%     %s1.A has all the contours and text position is   com(Aor(:,1:end),d1,d2);
%     c1 = plot_contours(s1.A,s1.V,0.95,1); %thr=.95,display=TRUE
%     fig1=gcf;
%      saveas(fig1,[figname1 '.fig'],'fig');
%      saveas(fig1,[figname1 '.png'],'png');
% end
kC1= fig1.Children.Children;
ncell=ceil((numel(kC1)-1)/2);
SIZ = size(fig1img);
SIZ = SIZ(1:2);
polygonfill=zeros(SIZ);
u=0;v=0;
F1Contour = nan(prod(SIZ),ncell);
Paninski = strcmp(class(kC1(1)),class(kC1(2)));
if Paninski
    for ii=1:ncell
        nth = ncell+1-ii;
        F1cell(ii).Position = [kC1(nth).Position(1)+v+5,kC1(nth).Position(2)+u+5];
        F1cell(ii).String = kC1(nth).String;
        F1Contour(:,ii) = reshape(circshift(kC1(nth+ncell).ZData,[u,v]),[],1);
    end
else
    for ii=1:ncell
        nth = 2*ncell+1-ii*2;
        F1cell(ii).Position = [kC1(nth).Position(1)+v+5,kC1(nth).Position(2)+u+5];
        F1cell(ii).String = kC1(nth).String;
        F1Contour(:,ii) = reshape(poly2mask(double(kC1(nth+1).XData),double(kC1(nth+1).YData),SIZ(1),SIZ(2)),[],1);
    end
end
    
%% select region for analysis
[f,p]=uigetfile('*.mat','Load region from previous file?');
if f~= 0
    load(fullfile(p,f),'maskpic0');
    %load(fullfile(p,f),'pair0');
else
    figure(fig1);
    maskpic0=roipoly;
end
close(fig1)
pair0 = reshape(maskpic0,1,prod(SIZ))*F1Contour >0;

%% load sig data -default for both eye 
% if menu('combined file for both eyes?','No','Yes') ==1
%     [eye1,eye2,sigtype]=loadsigdata; %sigtype: ca, sp, ML
%     SNR = cal_SNR;%%
% else
%     %eye1: ipsi, (stim file eye:2)
     [eye1,eye2,sigtype]=loadsigdata2; %sigtype: ca, sp, ML
     SNR = [];
% end    

plotrunrep(eye1.matrix,eye2.matrix)

%     Baseline substraction
    eye1.peak= eye1.peak-eye1.base;
    eye1.peakS=eye1.peakS-eye1.baseS;
    eye1.peakR=eye1.peakR-eye1.baseR;

    eye2.peak= eye2.peak-eye2.base;
    eye2.peakS=eye2.peakS-eye2.baseS;
    eye2.peakR=eye2.peakR-eye2.baseR;
    
for jj=1:size(eye1.peak,3)%the last one is background
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
    [a1,ori1]=nanmax(eye1.peak(:,:,nth1));%eye1:ipsi
    [a2,ori2]=nanmax(eye2.peak(:,:,nth2));%eye2:contra
    
    [a_s,pref_ori_s]=nanmax(nanmax(cat(1,eye1.peakS(:,:,nth1),eye2.peakS(:,:,nth2))));
    [a_r,pref_ori_r]=nanmax(nanmax(cat(1,eye1.peakR(:,:,nth1),eye2.peakR(:,:,nth2))));

    [a_1s,pref_ori_1s]=nanmax(eye1.peakS(:,:,nth1));
    [a_2s,pref_ori_2s]=nanmax(eye2.peakS(:,:,nth2));
    [a_1r,pref_ori_1r]=nanmax(eye1.peakR(:,:,nth1));
    [a_2r,pref_ori_2r]=nanmax(eye2.peakR(:,:,nth2));

    b_1s = eye1.baseS(:,pref_ori_1s,nth1);
    b_2s = eye2.baseS(:,pref_ori_2s,nth1);
    b_1r = eye1.baseR(:,pref_ori_1r,nth2);
    b_2r = eye2.baseR(:,pref_ori_2r,nth2);

    bsd_1s = eye1.baseSdS(:,pref_ori_1s,nth1);
    bsd_2s = eye2.baseSdS(:,pref_ori_2s,nth1);
    bsd_1r = eye1.baseSdR(:,pref_ori_1r,nth2);
    bsd_2r = eye2.baseSdR(:,pref_ori_2r,nth2);
    
    [odi,~] = ODIcalculation(pref_ori,eye2.peak(:,:,nth2),eye1.peak(:,:,nth1));
    [odi_s,~] = ODIcalculation(pref_ori,eye2.peakS(:,:,nth2),eye1.peakS(:,:,nth1));
    [odi_r,~] = ODIcalculation(pref_ori,eye2.peakR(:,:,nth2),eye1.peakR(:,:,nth1));
    [odi_ss,~] = ODIcalculation(pref_ori_s,eye2.peakS(:,:,nth2),eye1.peakS(:,:,nth1));
    [odi_rr,~] = ODIcalculation(pref_ori_r,eye2.peakR(:,:,nth2),eye1.peakR(:,:,nth1));

    ORI0(jj) = pref_ori;
    ORI1(:,jj) = [ori1 ori2];
    ORI2(:,jj) = [pref_ori_s pref_ori_1s pref_ori_2s pref_ori_r pref_ori_1r pref_ori_2r];
    ORI_1(:,jj) = [eye1.pref_ori eye2.pref_ori];
    ORI_2(:,jj) = [eye1.pref_ori_s eye2.pref_ori_s eye1.pref_ori_r eye2.pref_ori_r];
    
    PEAK0(jj) =a;%all responses
    PEAK00(:,jj)= [eye1.peak(:,pref_ori,nth1) eye1.peakS(:,pref_ori,nth1) eye1.peakR(:,pref_ori,nth1)...
        eye2.peak(:,pref_ori,nth2) eye2.peakS(:,pref_ori,nth2) eye2.peakR(:,pref_ori,nth2)]; % use combined pref_ori during run/still 
    PEAK1(:,jj)= [a1 a2];% eye-specific responses
    PEAK2(:,jj)= [a_s a_1s a_2s a_r a_1r a_2r]; %use pref_ori_s and pref_ori_r 
    BASELINE(:,jj) = [b_1s b_2s b_1r b_2r]; %use pref_ori_s and pref_ori_r 
    BASESD(:,jj)= [bsd_1s bsd_2s bsd_1r bsd_2r];
    
    ST0(:,jj)= [eye1.error(:,pref_ori,nth1)  eye2.error(:,pref_ori,nth2)];
    ST1(:,jj)= [eye1.errorS(:,pref_ori,nth1) eye2.errorS(:,pref_ori,nth2)...
        eye1.errorR(:,pref_ori,nth1) eye2.errorR(:,pref_ori,nth2)];
    ST2(:,jj)= [eye1.errorS(:,pref_ori_s,nth1) eye2.errorS(:,pref_ori_s,nth2)...
        eye1.errorR(:,pref_ori_r,nth1) eye2.errorR(:,pref_ori_r,nth2)];

    ODI0(jj) = odi;
    ODI1(:,jj) = [odi_s  odi_r]; %using same pref_ori
    ODI2(:,jj) = [odi_ss  odi_rr]; %using still/run pref_ori
    
    try
        RGB1(:,jj) = ODI2rgb(odi,a); %ODIcolorbar to make a reference barmap
    end
end 
%%
if ~exist('fign','var')
    fign = [fign0 sigtype '_'];
end
%%
% Set potential good cells criteria: response and size
try 
    pos_ = strfind(fign0,'_');
    sbxread(fullfile(path1a,fign0(1:pos_(3)-1)),0,1);
    global info;
    magnification = info.config.magnification;
    load('magnificationlist');
    mag= maglist(magnification)
catch
    disp('Update magnification otherwise mag=2'); 
    mag =2;
end
%%
c = calcontour(F1Contour,SIZ);
c.area_cutoff = mag^2*30 *[0.6 5] ; % cellsize = [96 396] under 2X
c.circ_cutoff = 1.6; % 1/3 < a/b < 3 c.circ=5/3
c.skew_cutoff = 1; 
  
pair1 = c.area<=c.area_cutoff(2) & c.area>=c.area_cutoff(1) ;
pair2 = c.circ<=c.circ_cutoff;
% responselimit = .03;
% pair1 = (sum(PEAK1)>=responselimit)&(SNR(1:end-1)>1.3)&pair1;
%
h21=figure();
imshow(fig1img,[])%ipsi eye,red;contra eye,green
hold on;
%goodcell red, bad
for jj=1:size(F1Contour,2)
    contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,jj),SIZ(1),SIZ(2)),[0.05 1],'LineColor',[pair0(jj) pair0(jj)&~pair1(jj) 0 ],'LineWidth',1);
%       text(F1cell(jj).Position(1),F1cell(jj).Position(2),[F1cell(jj).String sprintf('=%d,%.1f',c.area(jj),c.circ(jj))],...
%          'Color',[pair0(jj)&pair1(jj) 0 0],'FontSize',16)
end
try
    saveas(h21,[fign0 'cellshape.fig'])
end
saveas(h21,[fign0 'cellshape.png'])
%%
save([fign 'all.mat'],'SNR','eye1','eye2','ORI0','ORI1','ORI2','ORI_1','ORI_2',...
    'PEAK0','PEAK00','PEAK1','PEAK2','BASELINE','BASESD',...
    'ST0','ST1','ST2','ODI0','ODI1','ODI2',...
    'fig1img','F1cell','F1Contour','SIZ','RGB1','c',...
    'pair2','pair1','pair0','maskpic0','fign','fign0','mag','-v7.3');
disp(['saved: ' fign 'all.mat'])
%
selected1 = find(pair0&pair1);
% Analysis Summary for automatic selection
plotODIanalysis(fign,selected1)
%plotODIanalysis(fign,find(pair0&pair1& (sum(PEAK1)>prctile(sum(PEAK1),75))))
% plotODIanalysis(fign,find(pair0&~pair1))
close all
%%
%  plotTC(fign,selected1);
%  ManuelCellSelection
