% function analysis_twoday_data
% clear all;
close all;clc

[day1peakf,path1]=uigetfile('peakSI.mat','select peakSI matfile for day1eye');
[day2peakf,path2]=uigetfile('peakSI.mat','select peakSI matfile for day2eye');
% [alignedfig,figpath]=uigetfile('*_aligned.fig','find the aligned fig');
[alignedfig,figpath]=uigetfile('*&*.fig','find the aligned fig');
[prefilename,path]=uigetfile('*.mat','Pick the pre-matched pairs?');
%% load Day1 and Day2 traces and peakSI info
% [day1,day2,sigtype]=loadsigdata; %sigtype: ca, sp, ML


day1peak=load(fullfile(path1,day1peakf));
day1.finalvalue = permute(day1peak.SI.peakS,[3 2 1]);
day1.stdofeachvalue = permute(day1peak.SI.errorS,[3 2 1]);
day1.finalvalue2 = permute(day1peak.SI.peakR,[3 2 1]);
day1.stdofeachvalue2 = permute(day1peak.SI.errorR,[3 2 1]);
day2peak=load(fullfile(path2,day2peakf));
day2.finalvalue = permute(day2peak.SI.peakS,[3 2 1]);
day2.stdofeachvalue = permute(day2peak.SI.errorS,[3 2 1]);
day2.finalvalue2 = permute(day2peak.SI.peakR,[3 2 1]);
day2.stdofeachvalue2 = permute(day2peak.SI.errorR,[3 2 1]);

try
    day1.sigF= permute(day1peak.sigF,[4 3 2 1]);
    day2.sigF= permute(day2peak.sigF,[4 3 2 1]);
    suffix = '_CA';
catch
    day1.sigF= permute(day1peak.sigspF,[4 3 2 1]);
    day2.sigF= permute(day2peak.sigspF,[4 3 2 1]);
    suffix = '_SP';
end

% [day1sigf,path]=uigetfile('*signals.mat','select trace for day1eye');
% day1sig=load(fullfile(path,day1sigf))
% day1.sigF= permute(day1sig.sigF,[4 3 2 1]);
% [day2sigf,path]=uigetfile('*signals.mat','open trace for day2eye');
% day2sig=load(fullfile(path,day2sigf));
% day2.sigF= permute(day2sig.sigF,[4 3 2 1]);
%% sigF:ncells,nSteps, numberofcycles,framestocapture

nSteps=size(day1.sigF,2);
framestocapture=min (size(day1.sigF,4), size(day2.sigF,4));
nSteps1=nSteps-mod(nSteps,2);

if ~exist('F1Contour')
    pre=matfile(fullfile(path,prefilename));
    try
        pair1 = pre.pair1;
        pair2 = pre.pair2;
    catch
        pair1 = pre.selectedpair(:,1);
        pair2 = pre.selectedpair(:,2);
    end
    F1cell = pre.F1cell;
    F2cell = pre.F2cell;
    F1Contour = pre.F1Contour;
    F2Contour = pre.F2Contour;
end
selectedpair = cat(2,pair1(:),pair2(:));

%% Plot and pick cells with max_response
[names,names_dir] = num2ori(nSteps);
nCell = numel(pair1);
fit_selected = {};
peak_selected = nan(4,nSteps,nCell);
win_sig = day1peak.win_sig([1 end]);
%%
plotfig =1;
for j=1:numel(pair1)
    fit.peak1S = day1.finalvalue(pair1(j),:);
    [fit.a1S,fit.ori1S] = max(fit.peak1S);
    %     fit.ori1S = ori1s;
    % %     [fit.a1S,fit.ori1S,fit.peak1S(1:end-1)]=fit2peakgaussion(ori1s,fit.peak1S(1:end-1));    fit.peak1R = day1.finalvalue2(goodcell1(j),:);
    fit.peak1R = day1.finalvalue2(pair1(j),:);
    [fit.a1R,fit.ori1R] = max(fit.peak1R);
    %     [fit.a1R,fit.ori1R,fit.peak1R(1:end-1)]=fit2peakgaussion(ori1r,fit.peak1R(1:end-1));
    fit.peak2S = day2.finalvalue(pair2(j),:);
    [fit.a2S,fit.ori2S] = max(fit.peak2S);
    %     [fit.a2S,fit.ori2S,fit.peak2S(1:end-1)]=fit2peakgaussion(ori2s,fit.peak2S(1:end-1));
    %     plot(1:nSteps,fit.peak2S,'b-','Linewidth',2)
    fit.peak2R = day2.finalvalue2(pair2(j),:);
    [fit.a2R,fit.ori2R]= max(fit.peak2R);
    %     [fit.a2R,fit.ori2R,fit.peak2R(1:end-1)]=fit2peakgaussion(ori2r,fit.peak2R(1:end-1));
    fit_selected{j} = fit;
    peak_selected(:,:,j)=cat(1,fit.peak1S,fit.peak1R,fit.peak2S,fit.peak2R);
    temp_day1=squeeze(day1.sigF(pair1(j),:,:,:));
    temp_day2=squeeze(day2.sigF(pair2(j),:,:,:));
    yrange=prctile([temp_day1(:);temp_day2(:)],[0.5 99.5]);
    ymin = yrange(1);
    ymax = yrange(2);
    
    if plotfig 
    h1=figure('Position',[100 500 1000 450],'Name',['alltraces_cell#' num2str(pair1(j)) '_cell#' num2str(pair2(j))]);
    
    max_resp = nanmax([day1.finalvalue(pair1(j),:)+2*day1.stdofeachvalue(pair1(j),:),...
        day1.finalvalue2(pair1(j),:)+2*day1.stdofeachvalue2(pair1(j),:),...
        day2.finalvalue(pair2(j),:)+2*day2.stdofeachvalue(pair2(j),:),...
        day2.finalvalue2(pair2(j),:)+2*day2.stdofeachvalue(pair2(j),:)]);
    
try
    %plot peak amplitude
    subplot(2,3,1);    hold on;
     errorbar(1:nSteps,day1.finalvalue(pair1(j),:),day1.stdofeachvalue(pair1(j),:),'bo-');
    ori1s = day1peak.SI.pref_ori_S(pair1(j));
    %     plot(1:nSteps,fit.peak1S,'b--','Linewidth',2)
    errorbar(1:nSteps,day1.finalvalue2(pair1(j),:),day1.stdofeachvalue2(pair1(j),:),'ro-');
    ori1r = day1peak.SI.pref_ori_R(pair1(j));
    %     plot(1:nSteps,fit.peak1R,'r--','Linewidth',2)
    xlim([1 nSteps])
    ylim(yrange);
    title('Pre-exposure')  
    xlabel('orientation')
    ylabel('Response (df/F)');
    legend({'Still','Run'},'Orientation','horizontal')
    legend('boxoff')
    ax = gca;
    ax.XTickLabel = names_dir(ax.XTick);
%     legend([' Still ori=' num2str(ori1s)],['Run ori=' num2str(ori1r)])
    
     subplot(2,3,2);    hold on;
    errorbar(1:nSteps,day2.finalvalue(pair2(j),:),day2.stdofeachvalue(pair2(j),:),'o-','Color',[0 0 .5]);
    ori2s = day2peak.SI.pref_ori_S(pair2(j));
    
    errorbar(1:nSteps,day2.finalvalue2(pair2(j),:),day2.stdofeachvalue2(pair2(j),:),'o-','Color',[0.5 0 0]);
    ori2r = day2peak.SI.pref_ori_R(pair2(j));
    %     plot(1:nSteps,fit.peak2R,'r-','Linewidth',2)
    xlim([1 nSteps])
    ylim(yrange);
    %     legend(['Still ori=' num2str(ori2s)],[' ori=' num2str(fit.ori2S)],['Run ori=' num2str(ori2r)],[' ori=' num2str(fit.ori2R)] )
    title('Pose-exposure')  
    xlabel('orientation')
    ylabel('Response');
    legend({'Still','Run'},'Orientation','horizontal')
    legend('boxoff')
    ax = gca;
    ax.XTickLabel = names_dir(ax.XTick);
%     legend(['Still ori=' num2str(ori2s)],['Run ori=' num2str(ori2r)])
        
    subplot(2,3,3);    hold on;
    temp_peak= [fit.a1S fit.a1R fit.a2S fit.a2R];
    [~,pos] = max(temp_peak);
    temp_ori= [fit.ori1S fit.ori1R fit.ori2S fit.ori2R];
    fit.pref_ori = temp_ori(pos);
    plot(1:nSteps,fit.peak1S,'bo-','Linewidth',1)
    plot(1:nSteps,fit.peak1R,'ro-','Linewidth',1)
    plot(1:nSteps,fit.peak2S,'o-','Color',[ 0 0 .5],'Linewidth',1)
    plot(1:nSteps,fit.peak2R,'o-','Color',[ 0.5 0 0],'Linewidth',1);
%     title(['Pre&post exposure best ori=' num2str(fit.pref_ori)]);
    xlim([1 nSteps])
    ylim(yrange);
    title('Pre- vs. post- exposure');
    xlabel('orientation')
    ylabel('Response');
%    legend('Pre-','','Post-','' )
%     legend('boxoff')
    ax = gca;
    ax.XTickLabel = names_dir(ax.XTick);
end
   
    % Trace drawings

    
    for ori=1:nSteps
        traces1r = squeeze(temp_day1(ori,day1peak.matrix(:,ori),:));
        traces1s = squeeze(temp_day1(ori,~day1peak.matrix(:,ori),:));      
        traces2r = squeeze(temp_day2(ori,day2peak.matrix(:,ori),:));
        traces2s = squeeze(temp_day2(ori,~day2peak.matrix(:,ori),:));
      
        b=subplot(4,nSteps,2*nSteps+ori); hold on;
        rectangle('Position',[win_sig(1) 0 diff(win_sig) ymax],...
            'FaceColor',[0.5 .5 .5 .3],'EdgeColor',[0.5 .5 .5 .3]);% stimulus on&off time
%         if ori~=2 & ori~=4
            plot(traces1r','Color',[1 0 0 .3]);
            plot(squeeze(mean(temp_day1(ori,day1peak.matrix(:,ori),:))),'r','linewidth',2);
%         end
%         if ori~=1 & ori~=3
             plot(traces1s','Color',[0 0 1 .3]);
            plot(squeeze(mean(temp_day1(ori,~day1peak.matrix(:,ori),:))),'b','linewidth',2);
%         end
        xlim([1 framestocapture])
        ylim(yrange);
        
        %         title(num2str(day1.finalvalue2(goodcell1(j),ori)));
        Fr1(ori) = pwrplt(squeeze(mean(temp_day1(ori,:,:))),0);
%         title(sprintf('%.2f',Fr1(ori) ));
        
        c=subplot(4,nSteps,3*nSteps+ori); hold on;
        rectangle('Position',[win_sig(1) 0 diff(win_sig) ymax],...
            'FaceColor',[0.5 .5 .5 .3],'EdgeColor',[0.5 .5 .5 .3]);% stimulus on&off time
%         if ori~=2 & ori~=4
        plot(traces2r','Color',[.5 0 0 .3]);
        plot(squeeze(mean(temp_day2(ori,day2peak.matrix(:,ori),:))),'Color',[.5 0 0],'Linewidth',2);
%         end
%         if ori~=1 & ori~=3
            plot(traces2s','Color',[0 0 .5 .3]);
            plot(squeeze(mean(temp_day2(ori,~day2peak.matrix(:,ori),:))),'Color',[0 0 0.5],'Linewidth',2);
%         end
        xlim([1 framestocapture])
        ylim(yrange); 
        %         title(num2str(day2.finalvalue2(goodcell2(j),ori)))
        Fr2(ori) = pwrplt(squeeze(mean(temp_day2(ori,:,:))),0);
        title(names_dir(ori))
       if  ori ==ceil(nSteps/2)
            c.XLabel.String = 'Stimulus orientation';
       end
        if ori==1
            b.YLabel.String=['Pre-exposure']; %num2str(pair1(j)) ];
            c.YLabel.String=['Post-exposure'];  %num2str(pair2(j)) ];
            b.XTickLabel='';
            c.XTickLabel='';
        else
            b.YTickLabel='';
            c.YTickLabel='';
            b.XTickLabel='';
            c.XTickLabel='';
        end
    end
    
%     subplot(3,4,4);    hold on;
%     plot(1:nSteps,Fr1,'ko--','Linewidth',2)
%     plot(1:nSteps,Fr2,'ko-','Linewidth',2)
%     title(['FFT']);
%     xlim([1 nSteps])
%     legend('Baseline','Post-exposure' )
%     legend('boxoff')
   
    saveas(h1,[h1.Name '.fig'])
    saveas(h1,[h1.Name '.png'])
     close(h1)
    end
end
%%
filename=sprintf('selected_%dpairs_w_baseline',nCell);
save([filename '.mat'],'selectedpair','fit_selected','peak_selected');

%%
% hp = plotpeakselected(peak_selected);
% saveas(hp,[filename '.fig'])
% saveas(hp,[filename '.png'])

hp2 = plotpeakselected_env2(peak_selected);
%% Select cells
mkdir('Selected');
h4 = figure;
fig1imgnew = circshift(fig1img,[u,v]);
imshowpair(fig1imgnew,fig2img);
h4.Position = [1200 300 1200 1000];
SIZ = size(getimage(h4));

%%
selectedpair2 = zeros(size(selectedpair));
for j=1:nCell
    cell_1 = pair1(j);
    cell_2 = pair2(j);
    figure(h4);hold on
    ct1 = text(F1cell(cell_1).Position(1),F1cell(cell_1).Position(2),F1cell(cell_1).String,'Color',[ 1 0 0],'FontSize',8);
    cc1 = contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,cell_1),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[1 0 0]);
    ct2 = text(F2cell(cell_2).Position(1),F2cell(cell_2).Position(2),F2cell(cell_2).String,'Color',[ 0 1 0],'FontSize',8);
    cc2 = contour(1:SIZ(2),1:SIZ(1),reshape(F2Contour(:,cell_2),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[0 1 0]);
    
    figname = ['alltraces_cell#' num2str(cell_1) '_cell#' num2str(cell_2) '.png'];
    h5= figure('Position',[ 100 700 700 560]);
    A = imread(figname);
    imshow(A,'InitialMagnification',80)
    title(strrep(figname,'_','-'))
    if menu('Good cell?','Yes','No') ==1
            close(h5);
        selectedpair2(j,:)= selectedpair(j,:);
        copyfile(figname,'Selected');
    else
        close(h5);
        figure(h4);
        delete(ct1);
        delete(ct2);

    end
end
%%
fit_selected(selectedpair2(:,1)==0)=[];
peak_selected(:,:,selectedpair2(:,1)==0)=[];
selectedpair(selectedpair2(:,1)==0,:)=[];
newfolder = sprintf('selected_%dpairs',size(selectedpair,1));

filename= newfolder;
movefile('Selected',newfolder)
cd(newfolder)
save([filename '.mat'],'selectedpair','fit_selected','peak_selected');
saveas(h4,filename,'fig');

hp = plotpeakselected_env2(peak_selected);
saveas(hp,[filename '.fig'])
saveas(hp,[filename '.png'])

pair1 = selectedpair(:,1);
pair2 = selectedpair(:,2);

h0=figure('Name','preferred orientaion comparison');hold on;
preori1=[day1peak.SI.pref_ori_S(pair1);day1peak.SI.pref_ori_R(pair1)];
preori2=[day2peak.SI.pref_ori_S(pair2);day2peak.SI.pref_ori_R(pair2)];
arrow(preori1,preori2)
plot([0 max(preori1(:))],[0 max(preori2(:))],'--');
%legend('','day1','day2')
xlabel('Preferred orientation still');ylabel('Preferred orientation running');
saveas(h0,'preferred orientaion comparison.fig');

h1=figure('Name','gOSIcomparison');hold on;
gOSI1=squeeze([day1peak.SI.gOSI_S(pair1);day1peak.SI.gOSI_R(pair1)]);
gOSI2=squeeze([day2peak.SI.gOSI_S(pair2);day2peak.SI.gOSI_R(pair2)]);
plot([0 max(gOSI1(:))],[0 max(gOSI2(:))],'--');
xlabel('gOSI still');ylabel('gOSI running');
arrow(gOSI1,gOSI2)
saveas(h1,'gOSIcomparison.fig');

h2=figure('Name','OSIcomparison');hold on;
OSI1=squeeze([day1peak.SI.OSI_S(pair1);day1peak.SI.OSI_R(pair1)]);
OSI2=squeeze([day2peak.SI.OSI_S(pair2);day2peak.SI.OSI_R(pair2)]);
plot([0 max(OSI1(:))],[0 max(OSI2(:))],'--');
arrow(OSI1,OSI2)
xlabel('OSI still');ylabel('OSI running');
saveas(h2,'OSI running.fig');

h3=figure('Name','DSIcomparison');hold on;
DSI1=squeeze([day1peak.SI.DSI_S(pair1);day1peak.SI.DSI_R(pair1)]);
DSI2=squeeze([day2peak.SI.DSI_S(pair2);day2peak.SI.DSI_R(pair2)]);
plot([0 max(DSI1(:))],[0 max(DSI2(:))],'--');
arrow(DSI1,DSI2)
xlabel('DSI still');ylabel('DSI running');
saveas(h3,'DSIcomparison.fig');

%saveas(h,'OSI_twodays.fig')
%%
% day1=load('834_405_000&834_405_001_1_247cell_all');
% day2=load('834_410_001&834_410_002_1_263cell_all');
% for i=1:numel(pair1)
% cell1no(i)=str2num(F1cell(pair1(i)).String);
% cell2no(i)=str2num(F2cell(pair2(i)).String);
% end
% cell1no
% figure;title(selectedname)
% hold on;
% plot([day1.ODI1_bi(cell1no);day2.ODI1_bi(cell2no)])
% plot([mean(day1.ODI1_bi(cell1no));mean(day2.ODI1_bi(cell2no))],'ko-','Linewidth',2)
% xlim([.2 2.5])
