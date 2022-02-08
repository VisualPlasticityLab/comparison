
% function plot_twoday_data
% preselected pairs
%% load recorded ball,  cell data, and matched pairs
[fn1,folder1]=uigetfile('*ball.mat','select pre-stim ball matfile');
[fn2,folder2]=uigetfile('*ball.mat','select post-stim ball matfile');

load(fullfile(folder1,fn1),'velocity');vrun.pre=velocity(1:3:end);clear velocity
load(fullfile(folder2,fn2),'velocity');vrun.post=velocity(1:3:end);clear velocity
%% load cell data & matched pairs for all 3 planes
dd=dir('*selected*.mat');
for ii=1:3 
    
    allpairs{ii}=load(dd(ii).name);
    
    se1{ii}=load(fullfile(folder1,'suite2p',['plane' num2str(ii-1)],'Fall.mat'));
    se1{ii}.ops=RedChannelCorrection(se1{ii}.ops);  
    
    figpre1 = figure('Position',[10 300 100 600])
    plotRG(se1{ii},1:numel(se1{ii}.iscell));
    title(['plane' num2str(ii-1) 'pre-stim']);
    saveas(figpre1,['plane' num2str(ii-1) '.png'])
    
    figpre2 = figure('Position',[10 300 100 600])
    plotRG2(se1{ii},1:numel(se1{ii}.iscell));
    title(['plane' num2str(ii-1) 'pre-stim2']);
    saveas(figpre2,['plane' num2str(ii-1) '_2.png'])
    
    se2{ii}=load(fullfile(folder2,'suite2p',['plane' num2str(ii-1)],'Fall.mat'));
    se2{ii}.ops=RedChannelCorrection(se2{ii}.ops);
    
    figpost1 = figure('Position',[600 300 1000 600])
    plotRG(se2{ii},1:numel(se2{ii}.iscell));
    title(['plane' num2str(ii-1) 'post-stim']);
    saveas(figpost1,['plane' num2str(ii-1) '.png'])

    figpost2 = figure('Position',[600 300 1000 600])
    plotRG2(se2{ii},1:numel(se2{ii}.iscell));
    title(['plane' num2str(ii-1) 'post-stim2']);
    saveas(figpost2,['plane' num2str(ii-1) '_2.png'])
    
end
% %% calculate correlation 
% run_corr = corr(se1.velocity',se2.velocity');
% for n=1:numel(pair1)
%     sessions_corr(:,n) = corr(se1.F(pair1(n),:)',se2.F(pair2(n),:)');
%     run1_corr(:,n) = corr(se1.F(pair1(n),:)',se1.velocity');
%     run2_corr(:,n) = corr(se2.F(pair2(n),:)',se2.velocity');
% end
% 
% %% plot traces (modified from sigplt_onepage)
% 
% nframe=numel(se1.velocity);
% hsig=figure('Name','cell signals and correlation with running', 'position',[50 50 1500 50*n]);
% ax(1) = subplot(1,2,1)
% title('pre-stim')
% hold on;
% for n=1:numel(pair1)
%         Fsig = se1.F(pair1(n),:);
%         plot(1:nframe,Fsig/max(Fsig)*2+n,'LineWidth',1);   
%         text(nframe*1.00,n+.5,sprintf('#%d',pair1(n)),...
%             'HorizontalAlignment','left','VerticalAlignment','bottom');
% %           text(nframe*-.1,n,sprintf('C%d:%.2f',pair1(n),run1_corr(n)),...
% %             'HorizontalAlignment','right','VerticalAlignment','bottom');
% end
% plot(1:nframe,se1.velocity/max(se1.velocity)+n+2,'LineWidth',3,'Color',[0 0 0 .3]);
% text(nframe*1.00,n+2,'Run speed','HorizontalAlignment','left');
% axis off
% 
% ax(2) = subplot(1,2,2)
% title('post-stim')
% hold on;
% for n=1:numel(pair1)
%         Fsig = se2.F(pair2(n),:);
%         plot(1:nframe,Fsig/max(Fsig)*2+n,'LineWidth',1);   
%         text(nframe*1.00,n+.5,sprintf('#%d',pair2(n)),...
%             'HorizontalAlignment','left','VerticalAlignment','bottom');
% end
% plot(1:nframe,se2.velocity/max(se1.velocity)+n+2,'LineWidth',3,'Color',[0 0 0 .3]);
% text(nframe*1.00,n+2,'Run speed','HorizontalAlignment','left');
% axis off
% 
% linkaxes(ax)
% %% Plot correlation
% figure;hold on;
% plot([1;1.5;2], [run1_corr; sessions_corr;run2_corr],'x-','Color',[.5 .5 .5 .3],'LineWidth',.5);
% boxplot([run1_corr', run2_corr']);
% plot(1.5,sessions_corr)
% xlim([.5 2.5])
% %%
% filename=sprintf('selected_%dpairs',size(selectedpair,1));
% save([filename '.mat'],'selectedcell','fit_selected','peak_selected');
% 
% hp = plotpeakselected(peak_selected);
% saveas(hp,[filename '.fig'])
% saveas(hp,[filename '.png'])
% %%
% h0=figure('Name','preferred orientaion comparison');hold on;
% preori1=[se1.SI.pref_ori_S(goodcell1);se1.SI.pref_ori_R(goodcell1)];
% preori2=[se2.SI.pref_ori_S(goodcell2);se2.SI.pref_ori_R(goodcell2)];
% arrow(preori1,preori2)
% plot([0 max(preori1(:))],[0 max(preori2(:))],'--');
% %legend('','day1','day2')
% xlabel('Preferred orientation still');ylabel('Preferred orientation running');
% 
% h1=figure('Name','gOSIcomparison');hold on;
% gOSI1=squeeze([se1.SI.gOSI_S(goodcell1);se1.SI.gOSI_R(goodcell1)]);
% gOSI2=squeeze([se2.SI.gOSI_S(goodcell2);se2.SI.gOSI_R(goodcell2)]);
% plot([0 max(gOSI1(:))],[0 max(gOSI2(:))],'--');
% xlabel('gOSI still');ylabel('gOSI running');
% arrow(gOSI1,gOSI2)
% 
% h2=figure('Name','OSIcomparison');hold on;
% OSI1=squeeze([se1.SI.OSI_S(goodcell1);se1.SI.OSI_R(goodcell1)]);
% OSI2=squeeze([se2.SI.OSI_S(goodcell2);se2.SI.OSI_R(goodcell2)]);
% plot([0 max(OSI1(:))],[0 max(OSI2(:))],'--');
% arrow(OSI1,OSI2)
% xlabel('OSI still');ylabel('OSI running');
% 
% h3=figure('Name','DSIcomparison');hold on;
% DSI1=squeeze([se1.SI.DSI_S(goodcell1);se1.SI.DSI_R(goodcell1)]);
% DSI2=squeeze([se2.SI.DSI_S(goodcell2);se2.SI.DSI_R(goodcell2)]);
% plot([0 max(DSI1(:))],[0 max(DSI2(:))],'--');
% arrow(DSI1,DSI2)
% xlabel('DSI still');ylabel('DSI running');
% 
% %saveas(h,'OSI_twodays.fig')
% 
% %% Select cells
% mkdir('Selected');
% temp = figure;
% goodcell2 =  goodcell1;
% for j=nCell:-1:1
%     figname = ['alltraces_cell#' num2str(goodcell1(j)) '_cell#' num2str(goodcell2(j)) '.png'];
%     A = imread(figname);
%     imshow(A,'InitialMagnification',100),title(figname)
%     if menu('Good cell?','Yes','No') ==2
%         goodcell2(j) = [];
%         selectedpair(j,:)= 0;
% %         fit_selected(j) = [];
%         peak_selected(:,:,j)=[];
%     else
%         copyfile(figname,'Selected');
%     end
% end
% close(temp);
% 
% newfolder = sprintf('selected_%dpairs',size(selectedpair,1));
% movefile('Selected',newfolder)
% cd(newfolder)
% filename= newfolder;
% save([filename '.mat'],'selectedcell','fit_selected','peak_selected');
% hp = plotpeakselected_env(peak_selected);
% saveas(hp,[filename '.fig'])
% saveas(hp,[filename '.png'])
% 
% %%
% se1=load('834_405_000&834_405_001_1_247cell_all');
% se2=load('834_410_001&834_410_002_1_263cell_all');
% for n=1:numel(pair1)
% cell1no(n)=str2num(F1cell(pair1(n)).String);
% cell2no(n)=str2num(F2cell(pair2(n)).String);
% end
% cell1no
% figure;
% hold on;
% plot([se1.ODI1_bi(cell1no);se2.ODI1_bi(cell2no)])
% plot([mean(se1.ODI1_bi(cell1no));mean(se2.ODI1_bi(cell2no))],'ko-','Linewidth',2)
% xlim([.2 2.5])
% 
% %%# export them into a file, so all the individual neurons are being
% %%recorded and pulled together later to calculate the change for all mice

