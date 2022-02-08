%%%%%To be added: Fano factor for pre-screening?
%%%%%FFT selections?
clear all;close all;clc

%% load figures
[figname1,path1a]=uigetfile('.fig','find the fig and the memmap.mat');
fign=strtok(figname1,'t');
figimg = sbxshownplane;
if menu('Make sure 1st recording file is showed?','Yes','reselect') ==2
    figimg = sbxshownplane;
end
if iscell(figimg)
    fig1img = figimg{menu('which plane to select?','1st','2nd','3rd')};
else
    fig1img = figimg;
end

%% Select areas for cells

fig1=openfig(fullfile(path1a,figname1));
kC1= fig1.Children.Children;
ncell=(numel(kC1)-1)/2;
SIZ = size(fig1img);
SIZ = SIZ(1:2);
F1Contour =zeros(prod(SIZ),ncell);
pos = strfind(fign,'_');

u=0;v=0;
for ii=1:ncell
    nth = ncell+1-ii;
    F1cell(ii).Position = [kC1(nth).Position(1)+v+5,kC1(nth).Position(2)+u+5];
    F1cell(ii).String = kC1(nth).String;
    F1Contour(:,ii) = reshape(circshift(kC1(nth+ncell).ZData,[u,v]),[],1);
end
%% load data and fig
[~,ipsieyepath]=uigetfile('*.mat','select trace for ipsieye');
eye1=load(fullfile(ipsieyepath,'peakSI.mat'),'run','matrix','window','win_sig','sigF');
if isfield(eye1,'win_sig')
    win_sig1 = eye1.win_sig;
else
    load(fullfile(ipsieyepath,'peakSI.mat'),'window');
    win_sig1 = [eye1.window(1) eye1.window(end)];
end
if ~isfield(eye1,'matrix')
    load(fullfile(ipsieyepath,'peakSI.mat'),'run');
    eye1.matrix=run.matrix;
end
try
    [~,ipsieyepath2]=uigetfile('*.mat','Load more traces for ipsieye?');
    eye1extra=load(fullfile(ipsieyepath2,'peakSI.mat'),'sigF','matrix','win_sig');
    framestocapture=min(size(eye1extra.sigF,1),size(eye1.sigF,1));
    eye1.sigF = cat(2, eye1.sigF(1:framestocapture,:,:,:),eye1extra.sigF(1:framestocapture,:,:,:));
    if ~isfield(eye1extra,'matrix')
        load(fullfile(ipsieyepath2,'peakSI.mat'),'run');
        eye1extra.matrix=run.matrix;
    end
    eye1.matrix= cat(1,eye1.matrix,eye1extra.matrix);
end

[~,contraeyepath]=uigetfile('*.mat','select trace for contraeye');
eye2=load(fullfile(contraeyepath,'peakSI.mat'),'matrix','window','win_sig','sigF');
if isfield(eye2,'win_sig')
    win_sig2 = eye2.win_sig;
else
    load(fullfile(contraeyepath,'peakSI.mat'),'window');
    win_sig2 = [eye2.window(1) eye2.window(end)];
end
if ~isfield(eye2,'matrix')
    load(fullfile(contraeyepath,'peakSI.mat'),'run');
    eye2.matrix=run.matrix;
end
try
    [~,contraeyepath2]=uigetfile('*.mat','Load more traces for contraeye?');
    eye2extra=load(fullfile(contraeyepath2,'peakSI.mat'),'sigF','matrix');
    framestocapture=min(size(eye2extra.sigF,1),size(eye2.sigF,1));
    eye2.sigF = cat(2, eye2.sigF(1:framestocapture,:,:,:),eye2extra.sigF(1:framestocapture,:,:,:));
    if ~isfield(eye2extra,'matrix')
        load(fullfile(contraeyepath2,'peakSI.mat'),'run');
        eye2extra.matrix=run.matrix;
    end
    eye2.matrix= cat(1,eye2.matrix,eye2extra.matrix);
end

if isempty(eye1.matrix) | isempty(eye2.matrix)
    eye1.matrix = logical(ones(rep,nSteps));
    eye2.matrix = logical(ones(rep,nSteps));
end

[peakR1,sigR1,errorR1,peakS1,sigS1,errorS1]= sigFcmp(eye1.sigF,win_sig1,eye1.matrix);
[peakR2,sigR2,errorR2,peakS2,sigS2,errorS2]= sigFcmp(eye2.sigF,win_sig2,eye2.matrix);
[SI.peak,~,SI.error]=cal_ER(sigF,win_sig);  %peak:1*Var*ncell
        [SI.OSI,SI.pref_dir]=calOSI(SI.peak);
        SI.gOSI=calgOSI(SI.peak);
        SI.DSI=calDSI(SI.peak);
        
        
        %% Initialize parameters
pref_ori = nan(1,ncell);
ORI = nan(1,ncell);
ORI0= nan(1,ncell);
ORI1= zeros(1,ncell);
ODI = nan(1,ncell);
ODI1 = nan(1,ncell);
ODI0= nan(1,ncell);
ODI_bi = nan(1,ncell);
ODI1_bi = nan(1,ncell);
ODI0_bi= nan(1,ncell);
PEAK1= nan(1,ncell);
rep = size(eye1.sigF,2);
nSteps = size(eye1.sigF,3);
nSteps1 = nSteps-mod(nSteps,2);
framestocapture=min (size(eye1.sigF,1), size(eye2.sigF,1));
%% Set potential good cells criteria: response and size
responselimit = .05;
try 
    load(fullfile(path1a,[fign(1:pos(end-1)) 'memmap.mat']),'magnification');
catch
    magnification =2;
end
if mod(magnification,1)|| magnification==2
    magnification = magnification*2.5;
end 
load('magnificationlist');

cellsize = maglist(magnification)^2*30 *[1/3 3] ;
% cellsize = [50 450]
%% Calculate cells
h0=figure;
imshow(fig1img)
hold on;

for jj=1:ncell
    nth1= jj;
    nth2= jj;
    allpeaks = cat(1,peakR1(:,:,nth1),peakR2(:,:,nth2),peakS1(:,:,nth1),peakS2(:,:,nth2));
    [a1,ori1]=max(peakR1(:,:,nth1));
    [a2,ori2]=max(peakR2(:,:,nth2));
    [a,pref_ori]=nanmax(nanmax(cat(1,peakR1(:,:,nth1),peakR2(:,:,nth2))));
    [odi,odi_bi] = ODIcalculation(pref_ori,peakR2(:,:,nth2),peakR1(:,:,nth1));
    [~,pref_ori_s]=nanmax(nanmax(cat(1,peakS1(:,:,nth1),peakS2(:,:,nth2))));
    [odi_s,odi_bi_s] = ODIcalculation(pref_ori_s,peakS2(:,:,nth2),peakS1(:,:,nth1));
    ORI1(jj) = pref_ori+i*pref_ori_s;
    PEAK1(jj)= a;
    ODI1(jj) = odi;
    ODI1_s(jj) = odi_s;
    ODI1_bi(jj) = odi_bi;
    ODI1_bi_s(jj) = odi_bi_s;
    RGB1(:,jj) = ODI2rgb(odi,a); %ODIcolorbar to make a reference barmap
    % only draw potentially good cells within the selection region
end

pair1 = (PEAK1>=responselimit)& sum(F1Contour>.01)<=cellsize(2) & sum(F1Contour>.01)>=cellsize(1);
save([fign 'all.mat'],'ODI1','ODI1_bi','ODI1_s','ODI1_bi_s',...
    'fig1img','F1cell','RGB1','PEAK1','cellsize','pair1');


%% PLOTING individual TC
for jj=1:ncell
    nth1= jj;
    nth2= jj;
    if ( PEAK1(jj)>=responselimit) & (sum(F1Contour(:,nth1)>.01)<=cellsize(2)) & (sum(F1Contour(:,nth1)>.01)>=cellsize(1))
        ORI1(jj) = pref_ori+i*pref_ori_s;
        a = PEAK1(jj);
        odi=ODI1(jj) ;
        odi_s = ODI1_s(jj) ;
        odi_bi=ODI1_bi(jj);
        odi_bi_s=ODI1_bi_s(jj);
%             temp_1st=squeeze(eye1.sigspF(:,:,:,nth1));
%             temp_2nd=squeeze(eye2.sigspF(:,:,:,nth2));
        temp_1st=squeeze(eye1.sigF(:,:,:,nth1));
            temp_2nd=squeeze(eye2.sigF(:,:,:,nth2));

            ymax=prctile([temp_1st(:);temp_2nd(:)],99);
        ymin=0;
        buffer=(ymax-ymin)*.05;
        
        h1=figure('Position',[100 200 800 400],'Name',['alltraces_cell#' num2str(jj)]);
        for ori=1:nSteps
            traces1r = squeeze(temp_1st(:,eye1.matrix(:,ori),ori));
            traces1s = squeeze(temp_1st(:,~eye1.matrix(:,ori),ori));
            traces1r_avg=sigR1(:,ori,nth1);
            traces1s_avg=sigS1(:,ori,nth1);
            
            traces2r= squeeze(temp_2nd(:,eye2.matrix(:,ori),ori));
            traces2s= squeeze(temp_2nd(:,~eye2.matrix(:,ori),ori));
            traces2r_avg=sigR2(:,ori,nth2);
            traces2s_avg=sigS2(:,ori,nth2);
            
            try
                Fr1(ori)=pwrplt(squeeze(nanmean(temp_1st(:,:,ori),2)),1);
            catch
                Fr1(ori) = nan;
            end
            try
                Fr2(ori)=pwrplt(squeeze(nanmean(temp_2nd(:,:,ori),2)),1);
            catch
                Fr2(ori)= nan;
            end
            
            b=subplot(3,nSteps,nSteps+ori); hold on;
            plot(traces1r,'Color',[0 1 0 0.3]);
            plot(traces1r_avg,'Linewidth',2,'Color',[0 1 0 0.7]);
            
            plot(traces1s,'Color',[0 0 0 0.3]);
            plot(traces1s_avg,'Linewidth',1,'Color',[0 0 0 0.7]);
            plot([win_sig1(1) win_sig1(1)],[ymin ymax],'Linewidth',1,'Color',[0 0 0 0.3]);% stimulus on time
            plot([win_sig1(end) win_sig1(end)],[ymin ymax],'Linewidth',1,'Color',[0 0 0 0.3]);% stimulus off time
            title(sprintf('%.3f',peakR1(:,ori,nth1)));
            %         title(sprintf('%.2f',Fr1(ori))
            b.XLabel.String=sprintf('%d',ori);
            axis([1 framestocapture ymin-buffer ymax+buffer]);
            
            c=subplot(3,nSteps,2*nSteps+ori); hold on;
            plot(traces2s,'Color',[0 0 0 0.3]);
            plot(traces2s_avg,'Linewidth',1,'Color',[0 0 0 0.7]);
            
            plot(traces2r,'Color',[1 0 0 0.3]);
            plot(traces2r_avg,'Linewidth',2,'Color',[1 0 0 0.7]);
            plot([win_sig2(1) win_sig2(1)],[ymin ymax],'Linewidth',1,'Color',[0 0 0 0.3]);
            plot([win_sig2(end) win_sig2(end)],[ymin ymax],'Linewidth',1,'Color',[0 0 0 0.3]);
            axis([1 framestocapture ymin-buffer ymax+buffer]);
            title(sprintf('%.3f',peakR2(:,ori,nth2)));
            % title(sprintf('%.2f',Fr2(ori)))
            b.XTickLabel='';
            c.XTickLabel='';
            if ori==1
                b.YLabel.String=sprintf('ipsieye,cell# %d',nth1);
                c.YLabel.String=sprintf('contraeye,cell# %d',nth2);
            else
                b.YTickLabel='';
                c.YTickLabel='';
            end
        end
        
        subplot(3,3,1);    hold on;
        errorbar(1:nSteps,peakS1(:,:,nth1),errorS1(:,:,nth1),'Linewidth',1,'Color',[0 0 0 0.3]);
        errorbar(1:nSteps,peakR1(:,:,nth1),errorR1(:,:,nth1),'Linewidth',2,'Color',[0 1 0 0.3]);
        %         axis tight;
        title(['ipsi eye, ori=' num2str(ori1)] );
        
        subplot(3,3,2);    hold on;
        errorbar(1:nSteps,peakS2(:,:,nth2),errorS2(:,:,nth2),'Linewidth',1,'Color',[0 0 0 0.3]);
        errorbar(1:nSteps,peakR2(:,:,nth2),errorR2(:,:,nth2),'Linewidth',2,'Color',[1 0 0 0.3]);
        axis tight;
        title(['contra eye,ori=' num2str(ori2)] )
        
        subplot(3,3,3);    hold on;
        plot(1:nSteps,peakS1(:,:,nth1),'Linewidth',2,'Color',[0 1 0 0.3]);
        plot(1:nSteps,peakR1(:,:,nth1),'Linewidth',2,'Color',[0 1 0 0.7]);
        plot(1:nSteps,peakS2(:,:,nth2),'Linewidth',2,'Color',[1 0 0 0.3]);
        plot(1:nSteps,peakR2(:,:,nth2),'Linewidth',2,'Color',[1 0 0 0.7]);
        xlim([1 nSteps]);
        legend('ispi-R','ispi-S','contra-R','contra-S' );
        legend('ispi-R','contra-R');
        legend('boxoff','Location','top') ;
        title(sprintf('a=%.2f,ori=%d,odi=%.2f,odi_b_i=%.2f',a,pref_ori,odi,odi_bi));
        saveas(h1,[h1.Name '.png'])
        saveas(h1,[h1.Name '.fig'])
        close(h1)
    end
end
%% Select a ROI to print and analyze
[f,p]=uigetfile('*.mat','Load region from previous file?');
if f~= 0
    load(fullfile(p,f),'maskpic0');
    %load(fullfile(p,f),'pair0');
else
    figure(fig1);
    maskpic0=roipoly;
end

pair0 = reshape(maskpic0,1,prod(SIZ))*F1Contour >0;

%%
% selected = find(pair0&pair1&~out1);
selected = find(pair0&pair1);
h2=figure;
imshow(fig1img)
hold on;
%ipsi eye,red;contra eye,green
for ii =1:numel(ODI1)
    RGB1(:,ii) = ODI2rgb(ODI1(ii),PEAK1(ii));
end
for jj=selected
    contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,jj),SIZ(1),SIZ(2)),[0.01 1],'LineColor',RGB1(:,jj),'LineWidth',3)
end

ODI1_b=nanmean(ODI1_bi(selected));
ODI1_d=nanmean(ODI1(selected));
ODI1_b_s=nanmean(ODI1_bi_s(selected));
ODI1_d_s=nanmean(ODI1_s(selected));
numcell = num2str(numel(selected));

fprintf('overall %d cells, ODI= %.2f\n',ncell,nanmean(ODI1_bi));
fprintf('ROI overall %d cells, ODI= %.2f\n',sum(pair0),nanmean(ODI1_bi(pair0)));
fprintf('selected %s cells, ODI=%.2f\n',numcell,ODI1_b);

h3=figure;
seg=.3;
subplot(2,2,1);hold on;
histogram(ODI1_bi(selected),[-1.2:seg:1.2]);
title(sprintf('Total %s cell,ODI=%.2f',numcell,ODI1_b));

subplot(2,2,2);hold on;
histogram(ODI1(selected),[-1.2:seg:1.2]);
title(sprintf('ODI-direction =%.2f',ODI1_d));

subplot(2,2,3);hold on;
histogram(ODI1_bi_s(selected),[-1.2:seg:1.2]);
title(sprintf('ODI_still =%.2f',ODI1_b_s));

subplot(2,2,4);hold on;
histogram(ODI1_s(selected),[-1.2:seg:1.2]);
title(sprintf('ODI-direction_still =%.2f',ODI1_d_s));
%%
saveas(h2,[fign 'as' numcell '.fig']);
saveas(h2,[fign 'as' numcell '.png']);
saveas(h3,[fign 'as' numcell '_ODI.fig']);
saveas(h3,[fign 'as' numcell '_ODI.png']);
save([fign 'as' numcell],'maskpic0','pair0');

%% Manuel inspection
h22=figure;
imshow(fig1img)
hold on;
%ipsi eye,red;contra eye,green
pair2=zeros(1,ncell);
for jj=selected
     figure(h22);
     contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,jj),SIZ(1),SIZ(2)),[0.01 1],'LineColor',RGB1(:,jj),'LineWidth',3)
    imgtitle = ['alltraces_cell#' num2str(jj) '.png'];
    try
        img=imread(imgtitle);
        h1=figure;imshow(img);
        pair2(jj)= (menu('Good cell?','Y','N')==1);
        close(h1)
    catch
        pair2(jj) = 0;
    end
     contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,jj),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[1 1 1],'LineWidth',3)

end

selected = find(pair0&pair1&pair2);
h4=figure;
imshow(fig1img)%ipsi eye,red;contra eye,green
hold on;
for m=(selected)
    contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,jj),SIZ(1),SIZ(2)),[0.01 1],'LineColor',RGB1(:,jj),'LineWidth',3)
end
ODI1_b=nanmean(ODI1_bi(selected));
ODI1_d=nanmean(ODI1(selected));
ODI1_b_s=nanmean(ODI1_bi_s(selected));
ODI1_d_s=nanmean(ODI1_s(selected));
numcell = num2str(numel(selected));

h5=figure('Position',[ 100 200 600 800]);
subplot('Position',[0.1 0.1 0.8 0.6]);hold on
scatter(ODI1_s(selected),ODI1(selected),'o');
plot([-1 1],[-1 1],'--')
plot(mean(ODI1_s(selected)),mean(ODI1(selected)),'k*','MarkerSize',12);
ylabel('ODIrun')
xlabel('ODIstill')
% axis square

seg=.3;
h=subplot('Position',[0.1 0.8 0.8 .1]);hold on;
histplt2(ODI1_s(selected),ODI1(selected),[-1.2:seg:1.2]);
legend('ODIstill','ODIrun')
legend('Boxoff')
ylabel('Cell count')
title(['N=' numcell])
%%
for m=(selected)
    text(F1cell(jj).Position(1),F1cell(jj).Position(2),F1cell(jj).String,...
        'Color',RGB1(:,jj),'Color','k','FontSize',16)
end
redornot = F1Contour'*reshape(fig1img(:,:,1),[],1);
figure;histogram(redornot,0:20);
redornot(redornot<10)=0;

pair3 = redornot';
h4=figure;
imshow(fig1img)
hold on;
%ipsi eye,red;contra eye,green
for m=(selected)
    contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,jj),SIZ(1),SIZ(2)),[0.01 1],'LineColor',RGB1(:,jj),'LineWidth',3)
end
%%
for m=selected
    text(F1cell(jj).Position(1),F1cell(jj).Position(2),F1cell(jj).String,...
        'Color',RGB1(:,jj),'Color','k','FontSize',16)
end
%%
h4=figure; clf;
imshow(fig1img)
hold on;
%ipsi eye,red;contra eye,green
for m=selected
    contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,jj),SIZ(1),SIZ(2)),[0.01 1],...
        'LineColor',RGB1(:,jj),'Linewidth',3)
end
%%
ODI1_b=nanmean(ODI1_bi(selected));
ODI1_d=nanmean(ODI1(selected));
ODI1_b_s=nanmean(ODI1_bi_s(selected));
ODI1_d_s=nanmean(ODI1_s(selected));
numcell = num2str(numel(selected));

h5=figure('Position',[ 100 200 600 800]);
subplot('Position',[0.1 0.1 0.8 0.6]);hold on
scatter(ODI1_s(selected),ODI1(selected),'o');
plot([-1 1],[-1 1],'--')
plot(mean(ODI1_s(selected)),mean(ODI1(selected)),'k*','MarkerSize',12);
ylabel('ODIrun')
xlabel('ODIstill')
% axis square

seg=.3;
h=subplot('Position',[0.1 0.8 0.8 .1]);hold on;
histplt2(ODI1_s(selected),ODI1(selected),[-1.2:seg:1.2]);
legend('ODIstill','ODIrun')
legend('Boxoff')
ylabel('Cell count')
title(['N=' numcell])
% view([45 45 -90]);
% axis off

%%
save([fign 'as+ms' numcell],'maskpic0','pair0','pair2','pair3');

saveas(h4,[fign 'as+ms' numcell '.fig']);
saveas(h4,[fign 'as+ms' numcell '.png']);
saveas(h5,[fign 'as+ms' numcell '_ODI.fig']);
saveas(h5,[fign 'as+ms' numcell '_ODI.png']);
