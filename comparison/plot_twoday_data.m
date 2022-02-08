
% function plot_twoday_data
% preselected pairs
[day1,path1]=uigetfile('*_all.mat','select matfile for day1allcells');
[day2,path2]=uigetfile('*_all.mat','select matfile for day2allcells');
%% load Day1 and Day2 traces and peakSI info
% [day1sigf,path]=uigetfile('*signals.mat','select trace for day1eye');
% day1sig=load(fullfile(path,day1sigf))
% day1.sigF= permute(day1sig.sigF,[4 3 2 1]);

day1=load(fullfile(path1,day1));
day1.finalvalue = permute(day1.SI.peakS,[3 2 1]);
day1.stdofeachvalue = permute(day1.SI.errorS,[3 2 1]);
day1.finalvalue2 = permute(day1.SI.peakR,[3 2 1]);
day1.stdofeachvalue2 = permute(day1.SI.errorR,[3 2 1]);
day2=load(fullfile(path2,day2));
day2.finalvalue = permute(day2.SI.peakS,[3 2 1]);
day2.stdofeachvalue = permute(day2.SI.errorS,[3 2 1]);
day2.finalvalue2 = permute(day2.SI.peakR,[3 2 1]);
day2.stdofeachvalue2 = permute(day2.SI.errorR,[3 2 1]);

try
    day1.sigF= permute(day1.sigF,[4 3 2 1]);
    day2.sigF= permute(day2.sigF,[4 3 2 1]);
    suffix = '_CA';
catch
    day1.sigF= permute(day1.sigspF,[4 3 2 1]);
    day2.sigF= permute(day2.sigspF,[4 3 2 1]);
    suffix = '_SP';
end

%% sigF:ncells,nSteps, numberofcycles,framestocapture
nSteps=size(day1.sigF,2);
framestocapture=min (size(day1.sigF,4), size(day2.sigF,4));
nSteps1=nSteps-mod(nSteps,2);

[filename,path]=uigetfile('*.mat','Pick the pre-matched pairs?');
pre=load(fullfile(path,filename));
try
    goodcell1 = pre.pair1;
    goodcell2 = pre.pair2;
catch
    goodcell1 = pre.selectedcell(:,1);
    goodcell2 = pre.selectedcell(:,2);
end
matchedpair = cat(2,goodcell1,goodcell2);
%% Plot and pick cells with max_response
nCell = numel(goodcell1);
selectedpair= zeros(1,nCell);
fit_selected = {};
peak_selected = nan(4,nSteps,nCell);
%%
for j=1:1
    max_resp = cat(2,day1.finalvalue(goodcell1(j),:),day1.finalvalue2(goodcell1(j),:),day2.finalvalue(goodcell2(j),:),day2.finalvalue2(goodcell2(j),:));

    h1=figure('Position',[100 200 1200 900],'Name',['alltraces_cell#' num2str(goodcell1(j)) '_cell#' num2str(goodcell2(j))]);

    %plot peak amplitude
    subplot(3,4,1);    hold on;
    errorbar(1:nSteps,day1.finalvalue(goodcell1(j),:),day1.stdofeachvalue(goodcell1(j),:),'bo-');
    ori1s = day1.SI.pref_ori_S(goodcell1(j));
    fit.peak1S = day1.finalvalue(goodcell1(j),:);
    [fit.a1S,fit.ori1S] = max(fit.peak1S);
%     fit.ori1S = ori1s;
% %     [fit.a1S,fit.ori1S,fit.peak1S(1:end-1)]=fit2peakgaussion(ori1s,fit.peak1S(1:end-1));
%     plot(1:nSteps,fit.peak1S,'b--','Linewidth',2)
    errorbar(1:nSteps,day1.finalvalue2(goodcell1(j),:),day1.stdofeachvalue2(goodcell1(j),:),'ro-');
    ori1r = day1.SI.pref_ori_R(goodcell1(j));
    fit.peak1R = day1.finalvalue2(goodcell1(j),:);
    [fit.a1R,fit.ori1R] = max(fit.peak1R);
%     [fit.a1R,fit.ori1R,fit.peak1R(1:end-1)]=fit2peakgaussion(ori1r,fit.peak1R(1:end-1));
%     plot(1:nSteps,fit.peak1R,'r--','Linewidth',2)
    xlim([1 nSteps])
    title('Day1')
%     legend([' Still ori=' num2str(ori1s)],[' ori=' num2str(fit.ori1S)],['Run ori=' num2str(ori1r)] ,[' ori=' num2str(fit.ori1R)] )
    legend([' Still ori=' num2str(ori1s)],['Run ori=' num2str(ori1r)])

    legend('boxoff')
    
    subplot(3,4,2);    hold on;
    errorbar(1:nSteps,day2.finalvalue(goodcell2(j),:),day2.stdofeachvalue(goodcell2(j),:),'bo-');
    ori2s = day2.SI.pref_ori_S(goodcell2(j));
     fit.peak2S = day2.finalvalue(goodcell2(j),:);
     [fit.a2S,fit.ori2S] = max(fit.peak2S);
%     [fit.a2S,fit.ori2S,fit.peak2S(1:end-1)]=fit2peakgaussion(ori2s,fit.peak2S(1:end-1));
%     plot(1:nSteps,fit.peak2S,'b-','Linewidth',2)
    errorbar(1:nSteps,day2.finalvalue2(goodcell2(j),:),day2.stdofeachvalue2(goodcell2(j),:),'ro-');
    ori2r = day2.SI.pref_ori_R(goodcell2(j));
    fit.peak2R = day2.finalvalue2(goodcell2(j),:);
    [fit.a2R,fit.ori2R]= max(fit.peak2R);
%     [fit.a2R,fit.ori2R,fit.peak2R(1:end-1)]=fit2peakgaussion(ori2r,fit.peak2R(1:end-1));
%     plot(1:nSteps,fit.peak2R,'r-','Linewidth',2)
    xlim([1 nSteps])
    title('Day2')
%     legend(['Still ori=' num2str(ori2s)],[' ori=' num2str(fit.ori2S)],['Run ori=' num2str(ori2r)],[' ori=' num2str(fit.ori2R)] )
    legend(['Still ori=' num2str(ori2s)],['Run ori=' num2str(ori2r)])
    legend('boxoff')
    
    
    subplot(3,4,3);    hold on;
    temp_peak= [fit.a1S fit.a1R fit.a2S fit.a2R];
    [~,pos] = max(temp_peak);
    temp_ori= [fit.ori1S fit.ori1R fit.ori2S fit.ori2R];
    fit.pref_ori = temp_ori(pos);
    plot(1:nSteps,fit.peak1S,'bx--','Linewidth',2)
    plot(1:nSteps,fit.peak1R,'rx--','Linewidth',2)
    plot(1:nSteps,fit.peak2S,'bx-','Linewidth',2)
    plot(1:nSteps,fit.peak2R,'rx-','Linewidth',2); 
    title(['Day1&Day2 best ori=' num2str(fit.pref_ori)]);
    xlim([1 nSteps])
    legend('Day1 Still','Day1 Run','Day2 Still','Day2 run' )
    legend('boxoff')
    
% Trace drawings    
    temp_day1=squeeze(day1.sigF(goodcell1(j),:,:,:));
    temp_day2=squeeze(day2.sigF(goodcell2(j),:,:,:));
    ymax=max([temp_day1(:);temp_day2(:)]);
    ymin=min([temp_day1(:);temp_day2(:)]);
    buffer=(ymax-ymin)*.1;    
    drawnow;
    
        for ori=1:nSteps
        traces1 = squeeze(temp_day1(:,day1.matrix(:,ori),ori));
            traces1s = squeeze(temp_day2(:,~day1.matrix(:,ori),ori));
            traces1_avg=sigR1(:,ori,nth1);
            traces1s_avg=sigS1(:,ori,nth1);
            
            traces2= squeeze(temp_2nd(:,eye2.matrix(:,ori),ori));
            traces2s= squeeze(temp_2nd(:,~eye2.matrix(:,ori),ori));
            traces2_avg=sigR2(:,ori,nth2);
            traces2s_avg=sigS2(:,ori,nth2);
            
            
        b=subplot(3,nSteps,nSteps+ori); hold on;
        plot(squeeze(temp_day1(ori,:,:))');
        plot(squeeze(mean(temp_day1(ori,day1.matrix(:,ori),:))),'r','linewidth',3);
        plot(squeeze(mean(temp_day1(ori,~day1.matrix(:,ori),:))),'b','linewidth',3);
        plot([day1.win_sig(1) day1.win_sig(1)] ,[0 ymax] ,'k--','linewidth',1)
        plot([day1.win_sig(end) day1.win_sig(end)], [0 ymax],'k--','linewidth',1);% stimulus on time
        axis([1 framestocapture ymin-buffer ymax+buffer])
%         title(num2str(day1.finalvalue2(goodcell1(j),ori)));
        Fr1(ori) = pwrplt(squeeze(mean(temp_day1(ori,:,:))),0);
        title(sprintf('%.2f',Fr1));
        
        c=subplot(3,nSteps,2*nSteps+ori); hold on;
        plot(squeeze(temp_day2(ori,:,:))');
        plot(squeeze(mean(temp_day2(ori,day2.matrix(:,ori),:))),'r','linewidth',3);
        plot(squeeze(mean(temp_day2(ori,~day2.matrix(:,ori),:))),'b','linewidth',3);
        plot([day1.win_sig(1) day1.win_sig(1)] ,[0 ymax] ,'k--','linewidth',1)
        plot([day1.win_sig(end) day1.win_sig(end)], [0 ymax],'k--','linewidth',1);% stimulus on time
        axis([1 framestocapture ymin-buffer ymax+buffer])
%         title(num2str(day2.finalvalue2(goodcell2(j),ori)))
        Fr2(ori) = pwrplt(squeeze(mean(temp_day2(ori,:,:))),0);
        title(sprintf('%.2f',Fr2))
        if ori==1
            b.YLabel.String=['day1 #' num2str(goodcell1(j)) ];
            c.YLabel.String=['day2 #'  num2str(goodcell2(j)) ];
        else
            b.YTickLabel='';
            b.XTickLabel='';
            c.YTickLabel='';
            c.XTickLabel='';
        end
        end
    
    subplot(3,4,4);    hold on;
    plot(1:nSteps,Fr1,'ko--','Linewidth',2)
    plot(1:nSteps,Fr2,'ko-','Linewidth',2)
    title(['FFT']);
    xlim([1 nSteps])
    legend('Day1','Day2' )
    legend('boxoff')
    
%         fit_selected{j} = fit;
%         peak_selected(:,:,j)=cat(1,fit.peak1S,fit.peak1R,fit.peak2S,fit.peak2R);
%         saveas(h1,[h1.Name '.fig'])
%         saveas(h1,[h1.Name '.png'])
%         close(h1)
end
%%
filename=sprintf('selected_%dpairs',size(selectedpair,1));
save([filename '.mat'],'selectedcell','fit_selected','peak_selected');

%%
hp = plotpeakselected(peak_selected);
saveas(hp,[filename '.fig'])
saveas(hp,[filename '.png'])
%%
goodcell1 = selectedpair(:,1);
goodcell2 = selectedpair(:,2);
%%
h0=figure('Name','preferred orientaion comparison');hold on;
preori1=[day1.SI.pref_ori_S(goodcell1);day1.SI.pref_ori_R(goodcell1)];
preori2=[day2.SI.pref_ori_S(goodcell2);day2.SI.pref_ori_R(goodcell2)];
arrow(preori1,preori2)
plot([0 max(preori1(:))],[0 max(preori2(:))],'--');
%legend('','day1','day2')
xlabel('Preferred orientation still');ylabel('Preferred orientation running');

h1=figure('Name','gOSIcomparison');hold on;
gOSI1=squeeze([day1.SI.gOSI_S(goodcell1);day1.SI.gOSI_R(goodcell1)]);
gOSI2=squeeze([day2.SI.gOSI_S(goodcell2);day2.SI.gOSI_R(goodcell2)]);
plot([0 max(gOSI1(:))],[0 max(gOSI2(:))],'--');
xlabel('gOSI still');ylabel('gOSI running');
arrow(gOSI1,gOSI2)

h2=figure('Name','OSIcomparison');hold on;
OSI1=squeeze([day1.SI.OSI_S(goodcell1);day1.SI.OSI_R(goodcell1)]);
OSI2=squeeze([day2.SI.OSI_S(goodcell2);day2.SI.OSI_R(goodcell2)]);
plot([0 max(OSI1(:))],[0 max(OSI2(:))],'--');
arrow(OSI1,OSI2)
xlabel('OSI still');ylabel('OSI running');

h3=figure('Name','DSIcomparison');hold on;
DSI1=squeeze([day1.SI.DSI_S(goodcell1);day1.SI.DSI_R(goodcell1)]);
DSI2=squeeze([day2.SI.DSI_S(goodcell2);day2.SI.DSI_R(goodcell2)]);
plot([0 max(DSI1(:))],[0 max(DSI2(:))],'--');
arrow(DSI1,DSI2)
xlabel('DSI still');ylabel('DSI running');

%saveas(h,'OSI_twodays.fig')

%% Select cells
mkdir('Selected');
temp = figure;
goodcell2 =  goodcell1;
for j=nCell:-1:1
    figname = ['alltraces_cell#' num2str(goodcell1(j)) '_cell#' num2str(goodcell2(j)) '.png'];
    A = imread(figname);
    imshow(A,'InitialMagnification',100),title(figname)
    if menu('Good cell?','Yes','No') ==2
        goodcell2(j) = [];
        selectedpair(j,:)= 0;
%         fit_selected(j) = [];
        peak_selected(:,:,j)=[];
    else
        copyfile(figname,'Selected');
    end
end
close(temp);

newfolder = sprintf('selected_%dpairs',size(selectedpair,1));
movefile('Selected',newfolder)
cd(newfolder)
filename= newfolder;
save([filename '.mat'],'selectedcell','fit_selected','peak_selected');
hp = plotpeakselected_env(peak_selected);
saveas(hp,[filename '.fig'])
saveas(hp,[filename '.png'])

%%
day1=load('834_405_000&834_405_001_1_247cell_all');
day2=load('834_410_001&834_410_002_1_263cell_all');
for i=1:numel(pair1)
cell1no(i)=str2num(F1cell(pair1(i)).String);
cell2no(i)=str2num(F2cell(pair2(i)).String);
end
cell1no
figure;
hold on;
plot([day1.ODI1_bi(cell1no);day2.ODI1_bi(cell2no)])
plot([mean(day1.ODI1_bi(cell1no));mean(day2.ODI1_bi(cell2no))],'ko-','Linewidth',2)
xlim([.2 2.5])

%%# export them into a file, so all the individual neurons are being
%%recorded and pulled together later to calculate the change for all mice

