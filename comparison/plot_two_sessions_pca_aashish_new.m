
function plot_twoday_data(plane1, plane2, plane3)
% Aashish 
% preselected pairs
%% load recorded ball, cell data, and matched pairs
clear all
% [fn1,folder1]=uigetfile('*ball.mat','select pre-stim ball matfile');
% [fn2,folder2]=uigetfile('*ball.mat','select post-stim ball matfile');
% 
% load(fullfile(folder1,fn1),'velocity');pre.velocity=velocity(1:3:end);clear velocity
% load(fullfile(folder2,fn2),'velocity');post.velocity=velocity(1:3:end);clear velocity

%rather than manually selecting automatically select
%% Correct red channel based on green-red relationship
%alter dd based upon where the selected.mat files are found for each mouse
dd = dir('000-003\*selected*.mat');
% dd = dir('002-005\*selected*.mat');
split_file_name = strsplit(folder1, '\');
mouse_name = split_file_name{3};
exp_date = split_file_name{4};
for ii=1:3 
    % ii is number of the plane (1 to 3)
    
    allpairs{ii}=load(fullfile(dd(ii).folder, dd(ii).name));
    %contains matched cells 

    se{ii} = load(['F_' exp_date '_' mouse_name '_plane_' num2str(ii) '_JA.mat']);
    %loads file containing information about pre and post cells that are
    %matched across all 3 sessions
end
%% Generate redcell2 using Jennifer's 
for ii = 1:3
    se1{ii}=RedChannelCorrection(se1{ii});  
    se2{ii}=RedChannelCorrection(se2{ii});
end

%% plot all planes and label redcells based on orignial & JS's red
% for ii=1:3 
%    
%    f1 = plotRG(se1{ii},allpairs{ii}.pair1);    %default correction method
%     title(['plane' num2str(ii-1) 'pre-stim']);
%     
%     f2 =plotRG2(se1{ii},allpairs{ii}.pair1);
%     title(['plane' num2str(ii-1) 'pre-stim2']);
%         
%     f3 =plotRG(se2{ii},allpairs{ii}.pair2);
%     title(['plane' num2str(ii-1) 'post-stim']);
% 
%     f4= plotRG2(se2{ii},allpairs{ii}.pair2);
%     title(['plane' num2str(ii-1) 'post-stim2']);
%     
%     saveas(f1,['plane' num2str(ii-1) 'pre-stim'], 'png')
%     saveas(f2,['plane' num2str(ii-1) 'pre-stim2'], 'png')
%      saveas(f3,['plane' num2str(ii-1) 'post-stim'],'png')
%      saveas(f4,['plane' num2str(ii-1) 'post-stim2'],'png')
% 
% %   close all
% end
%% concatenate all data % matched pairs for all 3 planes - all 3 sessions 

%just f-fneu and stat
neur_coeff = 0.7;
pre.F =[];
pre.stat = []; 
post.F = [];
post.stat = [];

for ii=1:3
    pre.F = cat(1, pre.F, se{ii}.F_pre-neur_coeff*se{ii}.Fneu_pre);
    post.F = cat(1, post.F, se{ii}.F_post-neur_coeff*se{ii}.Fneu_post);
    pre.stat = cat(1,pre.stat,se{ii}.stat_pre);
    post.stat = cat(1,post.stat,se{ii}.stat_post);
end

npair = size(pre.F(1));
nframe = size(pre.F,2); 

%% concatenate all data & matched pairs for all 3 planes
% neur_coeff = 0.7;
% 
% pre.F =[];
% pre.spks = [];
% pre.redcell =[];
% pre.redcell2 =[];
% pre.stat = {};
% 
% for ii=1:3 
%     iscellidx = find(se1{ii}.iscell(:,1)==1)'; 
%     %find only the cells out of all ROI
%     pre.F = cat(1,pre.F, se1{ii}.F(iscellidx(allpairs{ii}.pair1),:)-neur_coeff*se1{ii}.Fneu(iscellidx(allpairs{ii}.pair1),:));
%     %getting the  calcium signal (neuropil corrected) only in the matched cells rather than all ROI
%     pre.spks = cat(1,pre.spks, se1{ii}.spks(iscellidx(allpairs{ii}.pair1),:));
%     %getting the deconvolved calcium signal (neuropil corrected) only in the matched cells rather than all ROI
%     pre.redcell= cat(1,pre.redcell, se1{ii}.redcell(iscellidx(allpairs{ii}.pair1),:));
%    pre.redcell2 = cat(1,pre.redcell2, se1{ii}.redcell2(iscellidx(allpairs{ii}.pair1),:));
%     pre.stat = cat(1,pre.stat,se1{ii}.stat{iscellidx(allpairs{ii}.pair1)});
% end
% 
% post.F =[];
% post.spks = [];
% post.redcell =[];
% post.redcell2 =[];
% post.stat = {};
% 
% for ii=1:3 
%     iscellidx = find(se2{ii}.iscell(:,1)==1)'; 
%     %find only the cells out of all ROI
%     post.F = cat(1,post.F, se2{ii}.F(iscellidx(allpairs{ii}.pair2),:)-neur_coeff*se2{ii}.Fneu(iscellidx(allpairs{ii}.pair2),:));
%     %getting the matched one from the cells rather than all ROI
%     post.spks = cat(1,post.spks, se2{ii}.spks(iscellidx(allpairs{ii}.pair2),:));
%     post.redcell= cat(1,post.redcell, se2{ii}.redcell(iscellidx(allpairs{ii}.pair2),:));
%    post.redcell2 = cat(1,post.redcell2, se2{ii}.redcell2(iscellidx(allpairs{ii}.pair2),:));
%     post.stat = cat(1,post.stat,se2{ii}.stat{iscellidx(allpairs{ii}.pair2)});
% end
% npair=size(pre.redcell,1);


%% plot traces (modified from sigplt_onepage)
% vmax = max(cat(2,pre.velocity,post.velocity));
% nframe=numel(pre.velocity);
% if npair>100
%     example= sort(unique(randi(npair,[1,100])));
% else 
%     example = 1:npair;
% end
% hsig=figure('Name','cell signals and correlation with running', 'position',[50 50 1500 50*numel(example)]);
% 
% 
% figaxes(1) = subplot(1,2,1);
% file_name = strsplit(se2{1}.ops.data_path, '/');
% mouse_name = file_name(3);
% % file_name = strsplit(se1{1}.ops.data_path, '/');
% % for i = 1:length(file_name)
% %     if file_name(i) == 'GAD*' 
% %         if strlength(file_name(i))== 4
% %             mouse_name = file_name(i)
% %         break
% %         end
% %     else 
% %         mouse_name = 'not found'
% %     end
% % end
% 
% %strfind(file_name, 'GAD')   
% 
% %title([mouse_name '  pre-stim'])
% title([mouse_name '  pre-stim']);
% hold on;
% ii=0;
% for n=example
%     ii=ii+1;
%         Fsig = pre.F(n,:);
%         plot(1:nframe,Fsig/max(Fsig)*2+ii,'LineWidth',1);   
%         text(nframe*1.00,ii+.5,sprintf('#%d',n),...
%             'HorizontalAlignment','left','VerticalAlignment','bottom');
% %           text(nframe*-.1,n,sprintf('C%d:%.2f',pair1(n),run1_corr(n)),...
% %             'HorizontalAlignment','right','VerticalAlignment','bottom');
% end
% plot(1:nframe,pre.velocity+ii+2,'LineWidth',3,'Color',[0 0 0 .3]);
% text(nframe*1.00,ii+2,'Run speed','HorizontalAlignment','left');
% axis off
% line([1 1],[ii+2 ii+8],'Linewidth',3)
% %strfind, 
% figaxes(2) = subplot(1,2,2);
% %title([mouse_name '  pre-stim'])
% title([mouse_name '  post-stim']);
% hold on;
% ii=0;
% for n=example
%             ii=ii+1;
%         Fsig = post.F(n,:);
%         plot(1:nframe,Fsig/max(Fsig)*2+ii,'LineWidth',1);   
%         text(nframe*1.00,ii+.5,sprintf('#%d',n),...
%             'HorizontalAlignment','left','VerticalAlignment','bottom');
% end
% plot(1:nframe,post.velocity+ii+2,'LineWidth',3,'Color',[0 0 0 .3]);
% text(nframe*1.00,ii+2,'Run speed','HorizontalAlignment','left');
% axis off
% line([1 1],[ii+2 ii+8],'Linewidth',3)
% 
% linkaxes(figaxes);

%% calculate & plot correlation
% 
% run_corr = corr(pre.velocity',post.velocity');
% %user set threshold - function - input threshold or default max velocity
% %from pre and post
% max_velocity_both(1) = max(pre.velocity);
% max_velocity_both(2) = max(post.velocity);
% max_velocity = sort(max_velocity_both, 'descend');
% prompt = 'Set a threshold for velocity in cm/s. If not necessary press ENTER and default value of maximum velocity will be assigned: ';
% thr = input(prompt);
% %try catch
% %ONE SET
% while isnumeric(thr) == 0
%     if isempty(thr);
%         break
%     else
%         prompt = 'Set a threshold for velocity in cm/s. If not necessary press ENTER and default value of maximum velocity will be assigned: ';
%         thr = input(prompt);
%     end  
% end
% if isempty(thr);
%     thr = max_velocity(1) + 1;
% end
% if thr < 1
%     thr_status = 'still';
% else
%     thr_status = 'full';
% end
% %ONE SET
% % if str == 'N'
% %     thr = max_velocity(1)
% % elseif str == 'Y'
% %     thr = str2double(input('Please enter your threshold in cm/s: '))
% % elseif isempty(str)
% %     thr = max_velocity(1)    
% % else
% %     thr = max_velocity(1)
% % end
% 
% %thr=0.05; % max running speed to be compared, if 0.05 then similar to still
% % nonrun=find(pre.velocity<=thr&post.velocity<=thr);
% 
% section=(pre.velocity<=thr&post.velocity<=thr);
% sect = (diff(section)==1);
% % nonrun_sect = 
% for n=1:npair
%     run1_corr(n,:) = corr(pre.F(n,section)',pre.velocity(section)');
%     run2_corr(n,:) = corr(post.F(n,section)',post.velocity(section)');
%     sessions_corr(n,:) = corr(pre.F(n,section)',post.F(n,section)');
% end
% % Plot correlation
% figure;
% boxplot([run1_corr, sessions_corr, run2_corr]);
% hold on;
% plot([1;2;3], [run1_corr, sessions_corr,run2_corr],'x-','Color',[.5 .5 .5 .3],'LineWidth',.5);
% %plot(1.5,sessions_corr)
% xlim([.5 3.5])
% xticks(1:3)
% xticklabels({'run-pre','pre-post','run-post'})
% title(sprintf('correlation-allpairs,runcorr=%.2f',run_corr))
% 
%% plot skewness
% for i = 1:numel(pre.stat);
%     pre_skew(i) = pre.stat{i}.skew;
% end
% for j = 1:numel(post.stat);
%     post_skew(j) = post.stat{j}.skew;
% end
% 
% figure; hold on;
% boxplot([pre_skew', post_skew']);
% plot([1;2], [pre_skew', post_skew'],'x-', 'Color',[.5 .5 .5 .3], 'LineWidth',.5);
% title([pre.ops.mouse_name ' skewness, pre and post'])
% xlim([.5 2.5]);
% xticks(1:2);
% xticklabels({'pre-stimulation-skew','post-stimulation-skew'})
%% Edit movie frames
movie_frames_pre = [];
movie_frames_post = [];

if se{1}.ops_pre.start_frame >= se{1}.ops_post.start_frame
    start_frame_later = se{1}.ops_pre.start_frame;
else
    start_frame_later = se{1}.ops_post.start_frame;
end

movie_sect_length = nframe - start_frame_later + 1;
start_frame_pre = se{1}.ops_pre.start_frame;
start_frame_post = se{1}.ops_post.start_frame;
end_frame_pre = se{1}.ops_pre.start_frame + movie_sect_length - 1;
end_frame_post = se{1}.ops_post.start_frame + movie_sect_length - 1;

for i = 1:size(pre.F,2)
    if i >= start_frame_pre & i<= end_frame_pre
        movie_frames_pre(end+1) = i;
    end
end
for i = 1:size(pre.F,2)
    if i >= start_frame_post & i<= end_frame_post
        movie_frames_post(end+1) = i;
    end
end
    
%% pca analysis
%z-score the data
%run pca analysis and plot
%do it on prestim and see how it correlates with running
%look at poststim data with the same set of dimensions from pca on session1
%use above code to plot

%normalising the data

normalised_pre = zscore(pre.F(:,movie_frames_pre), 0, 2);
normalised_post = zscore(post.F(:,movie_frames_post), 0, 2);
normalised = cat(2,normalised_pre,normalised_post);
%run pca analysis

[coeff, score, latent, ~, explained] = pca(normalised');
% [coeff_1, score_1, latent_1, ~, explained_1] = pca(normalised_pre');
% [coeff_2, score_2, latent_2, ~, explained_2] = pca(normalised_post');

%eigenvectors, scores(remapping of data), eigenvalues, ~, percentage of
%variance explained by the PC

% %only for my PCA, it has been done in matlab pca function
% normalised_pre0=normalised_pre; %-mean(normalised_pre,2);
% normalised_post0=normalised_post; %-mean(normalised_post,2);
% normalised0 = normalised; % - mean(normalised,2);


c = turbo(movie_sect_length);

%BOTH SESSIONS
figure;hold on;
title('PCA scores');
%scores
subplot(1,2,1); hold on;
for i = 1:movie_sect_length
    plot3(score(i,1), score(i,2), score(i,3), '.', 'color', c(i,:))

end
figaxes(1) = gca;
figaxes(1).XLabel.String = 'PC1';
figaxes(1).YLabel.String = 'PC2';
figaxes(1).ZLabel.String = 'PC3';
figaxes(1).Title.String = [mouse_name ' pre-stim'] % thr_status];
view(135, 20);
grid on;
subplot(1,2,2); hold on;
for i = (movie_sect_length + 1):(movie_sect_length + movie_sect_length)
    plot3(score(i,1), score(i,2), score(i,3), '.', 'color', c(i - movie_sect_length,:))
    
end
figaxes(2) = gca;
figaxes(2).XLabel.String = 'PC1';
figaxes(2).YLabel.String = 'PC2';
figaxes(2).ZLabel.String = 'PC3';
figaxes(2).Title.String = [mouse_name ' post-stim'] % thr_status];
view(135, 20);
axes_array_1(1) = figaxes(1).ZLim(1);
axes_array_1(2) = figaxes(2).ZLim(1);
axes_array_2(1) = figaxes(1).ZLim(2);
axes_array_2(2) = figaxes(2).ZLim(2);
axes_array_sorted_1 = sort(axes_array_1);
axes_array_sorted_2 = sort(axes_array_2);
Z_Lim(1) = axes_array_sorted_1(1);
Z_Lim(2) = axes_array_sorted_2(2);
figaxes(1).ZLim = Z_Lim;
figaxes(2).ZLim = Z_Lim;
grid on;

linkaxes(figaxes);
color_bar = colorbar;
color_bar.Ticks = [0 1];
color_bar.TickLabels = {'1', num2str(movie_sect_length)};
color_bar.Label.String = 'Movie Frame Num';
colormap turbo;

%% calculate correlation

for i = 1:size(pre.F,1)
    R = corrcoef(score(1:movie_sect_length,i), score((movie_sect_length + 1):movie_sect_length * 2,i));
    pca_corr_coeff(i) = R(1,2);
end

figure;hold on;
for i = 1:size(pre.F,1)
    plot(i, pca_corr_coeff(i), 'x', 'color', 'r')
end
title([mouse_name ': corr pre vs post for each PC']);
xlabel('PC');
ylabel('correlation_coefficient');

%%
% %vectors
% for i = 1:npair;
%     plot([0, Vec_sorted(i,1)],[0, Vec_sorted(i,2)], 'color', 'r')
%     text(Vec_sorted(i,1), Vec_sorted(i,2), num2str(i))
% end


% %ONE SESSION - USING MY GENERATION OF SCORES
% %prestim mine
% figure;hold on;
% title('PCA scores, prestim');
% %scores
% for i = 1:nframe;
%     plot(cov_test_score_pre(i,1), cov_test_score_pre(i,2), '.', 'color', c(i,:))
% end
% for i = nframe+1:nframe*2
%     plot(cov_test_score_pre(i,1), cov_test_score_pre(i,2), 'x', 'color', c(i,:)/2)
% end
% color_bar = colorbar;
% % color_bar.Tick
% % Labels=(0 1550)
% colormap turbo;
% 
% 
% 
% %poststim mine
% figure;hold on;
% title('PCA scores, poststim');
% for i = 1:nframe;
%     plot(cov_test_score_post(i,1), cov_test_score_post(i,2), '.', 'color', c(i,:))
% end
% for i = 1:npair;
%     plot([0, V_sorted(i,1)],[0, V_sorted(i,2)], 'color', 'b')
%     text(V_sorted(i,1), V_sorted(i,2), num2str(i))
% end
% color_bar = colorbar;
% % color_bar.Tick
% % Labels=(0 1550)
% colormap turbo
% 
% %NOT USING MY GENERATION OF SCORES
% %prestim
% figure;hold on;
% title('PCA scores, prestim');
% 
% %scores
% for i = 1:nframe;
%     %biplot(coeff(:,1:2),'scores',score(:,1:2),)
%     plot(score_1(i,1), score_1(i,2), '.', 'color', c(i,:))
% end
% 
% color_bar = colorbar;
% color_bar.Limits = [0 nframe];
% 
% %vectors
% for i = 1:npair;
%     plot([0, coeff_1(i,1)],[0, coeff_1(i,2)], 'color', 'b')
%     text(coeff_1(i,1), coeff_1(i,2), num2str(i))
% end
% 
% %colour
% figure; hold on;
% for cc = 1:nframe;
%     plot(cc, cc, '.', 'color', c(cc,:));
% end
% 
% %plot with no movement





%for the biplot one
% for i = 1:numel(pre.F(:,1))    
%     a{i} = num2str(i)
% end
% c = turbo(nframe);
% figure;hold on;
% for i = 1:nframe;
% %    biplot(coeff_1(:,1:2),'scores',score_1(:,1:2));
% end
%,'varlabels',{a'});

%%
% filename=sprintf('selected_%dpairs',npair);
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
% 
