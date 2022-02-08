
% function plot_twoday_data
% preselected pairs
%% load recorded ball, cell data, and matched pairs
[fn1,folder1]=uigetfile('*ball.mat','select pre-stim ball matfile');
[fn2,folder2]=uigetfile('*ball.mat','select post-stim ball matfile');

load(fullfile(folder1,fn1),'velocity');pre.velocity=velocity(1:3:end);clear velocity
load(fullfile(folder2,fn2),'velocity');post.velocity=velocity(1:3:end);clear velocity
%% Correct red channel based on green-red relationship
dd=dir('*selected*.mat');

for ii=1:3 
    % ii is number of the plane (1 to 3)
    
    allpairs{ii}=load(dd(ii).name);
    %contains matched cells 
    se1{ii}=load(fullfile(folder1,'suite2p',['plane' num2str(ii-1)],'Fall.mat'));
    %se1 is for prestim 
    se2{ii}=load(fullfile(folder2,'suite2p',['plane' num2str(ii-1)],'Fall.mat'));
    %se2 is for poststim
    
end
%%
for ii = 1:3
    se1{ii}=RedChannelCorrection(se1{ii});  
    se2{ii}=RedChannelCorrection(se2{ii});
end

%% plot all planes and label redcells based on orignial & JS's red
for ii=1:3 
   
   f1 = plotRG(se1{ii},allpairs{ii}.pair1);    %defult correction method
    title(['plane' num2str(ii-1) 'pre-stim']);
    
    f2 =plotRG2(se1{ii},allpairs{ii}.pair1);
    title(['plane' num2str(ii-1) 'pre-stim2']);
        
    f3 =plotRG(se2{ii},allpairs{ii}.pair2);
    title(['plane' num2str(ii-1) 'post-stim']);

    f4= plotRG2(se2{ii},allpairs{ii}.pair2);
    title(['plane' num2str(ii-1) 'post-stim2']);
    
    saveas(f1,['plane' num2str(ii-1) 'pre-stim'], 'png')
    saveas(f2,['plane' num2str(ii-1) 'pre-stim2'], 'png')
     saveas(f3,['plane' num2str(ii-1) 'post-stim'],'png')
     saveas(f4,['plane' num2str(ii-1) 'post-stim2'],'png')

%   close all
end

%% concatenate all data & matched pairs for all 3 planes
    
pre.F =[];
pre.spks = [];
pre.redcell =[];
pre.redcell2 =[];
pre.stat = {};

for ii=1:3 
    iscellidx = find(se1{ii}.iscell(:,1)==1)'; 
    %find only the cells out of all ROI
    pre.F = cat(1,pre.F, se1{ii}.F(iscellidx(allpairs{ii}.pair1),:)-se1{ii}.ops.neucoeff*se1{ii}.Fneu(iscellidx(allpairs{ii}.pair1),:));
    %getting the  calcium signal (neuropil corrected) only in the matched cells rather than all ROI
    pre.spks = cat(1,pre.spks, se1{ii}.spks(iscellidx(allpairs{ii}.pair1),:));
    %getting the deconvolved calcium signal (neuropil corrected) only in the matched cells rather than all ROI
    pre.redcell= cat(1,pre.redcell, se1{ii}.redcell(iscellidx(allpairs{ii}.pair1),:));
%    pre.redcell2 = cat(1,pre.redcell2, se1{ii}.redcell2(iscellidx(allpairs{ii}.pair1),:));
    pre.stat = cat(1,pre.stat,se1{ii}.stat{iscellidx(allpairs{ii}.pair1)});
end

post.F =[];
post.spks = [];
post.redcell =[];
post.redcell2 =[];
post.stat = {};

for ii=1:3 
    iscellidx = find(se2{ii}.iscell(:,1)==1)'; 
    %find only the cells out of all ROI
    post.F = cat(1,post.F, se2{ii}.F(iscellidx(allpairs{ii}.pair2),:)-se2{ii}.ops.neucoeff*se2{ii}.Fneu(iscellidx(allpairs{ii}.pair2),:));
    %getting the matched one from the cells rather than all ROI
    post.spks = cat(1,post.spks, se2{ii}.spks(iscellidx(allpairs{ii}.pair2),:));
    post.redcell= cat(1,post.redcell, se2{ii}.redcell(iscellidx(allpairs{ii}.pair2),:));
%    post.redcell2 = cat(1,post.redcell2, se2{ii}.redcell2(iscellidx(allpairs{ii}.pair2),:));
    post.stat = cat(1,post.stat,se2{ii}.stat{iscellidx(allpairs{ii}.pair2)});
end
npair=size(pre.redcell,1);


%% plot traces (modified from sigplt_onepage)
vmax = max(cat(2,pre.velocity,post.velocity));
nframe=numel(pre.velocity);
if npair>100
    example= sort(unique(randi(npair,[1,100])));
else 
    example = 1:npair;
end
hsig=figure('Name','cell signals and correlation with running', 'position',[50 50 1500 50*numel(example)]);


ax(1) = subplot(1,2,1)

title([fn1(1:4) '  pre-stim'])
hold on;
ii=0;
for n=example
    ii=ii+1;
        Fsig = pre.F(n,:);
        plot(1:nframe,Fsig/max(Fsig)*2+ii,'LineWidth',1);   
        text(nframe*1.00,ii+.5,sprintf('#%d',n),...
            'HorizontalAlignment','left','VerticalAlignment','bottom');
%           text(nframe*-.1,n,sprintf('C%d:%.2f',pair1(n),run1_corr(n)),...
%             'HorizontalAlignment','right','VerticalAlignment','bottom');
end
plot(1:nframe,pre.velocity+ii+2,'LineWidth',3,'Color',[0 0 0 .3]);
text(nframe*1.00,ii+2,'Run speed','HorizontalAlignment','left');
axis off
line([1 1],[ii+2 ii+8],'Linewidth',3)
%strfind, 
ax(2) = subplot(1,2,2)
title([fn1(1:4) '  post-stim'])
hold on;
ii=0;
for n=example
            ii=ii+1;
        Fsig = post.F(n,:);
        plot(1:nframe,Fsig/max(Fsig)*2+ii,'LineWidth',1);   
        text(nframe*1.00,ii+.5,sprintf('#%d',n),...
            'HorizontalAlignment','left','VerticalAlignment','bottom');
end
plot(1:nframe,post.velocity+ii+2,'LineWidth',3,'Color',[0 0 0 .3]);
text(nframe*1.00,ii+2,'Run speed','HorizontalAlignment','left');
axis off
line([1 1],[ii+2 ii+8],'Linewidth',3)

linkaxes(ax)

%% calculate & plot correlation

run_corr = corr(pre.velocity',post.velocity');
thr=0.05; % max running speed to be compared, if 0.05 then similar to still
% nonrun=find(pre.velocity<=thr&post.velocity<=thr);
section=(pre.velocity<=thr&post.velocity<=thr);
sect = (diff(section)==1);
% nonrun_sect = 
for n=1:npair
    run1_corr(n,:) = corr(pre.spks(n,section)',pre.velocity(section)');
    run2_corr(n,:) = corr(post.spks(n,section)',post.velocity(section)');
    sessions_corr(n,:) = corr(pre.spks(n,section)',post.spks(n,section)');
end
% Plot correlation
figure;
boxplot([run1_corr, sessions_corr, run2_corr]);
hold on;
plot([1;2;3], [run1_corr, sessions_corr,run2_corr],'x-','Color',[.5 .5 .5 .3],'LineWidth',.5);
%plot(1.5,sessions_corr)
xlim([.5 3.5])
xticks(1:3)
xticklabels({'run-pre','pre-post','run-post'})
title(sprintf('correlation-allpairs,runcorr=%.2f',run_corr))

%% pca analysis
%z-score the data
%run pca analysis and plot
%do it on prestim and see how it correlates with running
%look at poststim data with the same set of dimensions from pca on session1
%use above code to plot

%normalising the data

normalised_pre = zscore(pre.F, 0, 1);
normalised_post = zscore(post.F, 0, 1);

pre_minus_mean = [];
for iii = 1:npair
    row_mean(iii) = (sum(pre.F(iii,:)))/numel(pre.F(iii,:));
    pre_minus_mean(iii, :) = pre.F(iii,:) - row_mean(iii);
end


%run pca analysis

[coeff, score, latent] = pca(testing_data')
[coeff_1, score_1, latent_1, ~, explained_1] = pca(normalised_pre');
[coeff_2, score_2, latent_2, ~, explained_2] = pca(normalised_post');

loadings_1 = V_sorted(:,1) * sqrt(var(score_1(:,1)));
loadings_coeff_1 = coeff_1(:,1) 
test_score = loadings_1 * normalised_pre;
%covariance matrix, eigenvectors, eigenvalues
cov_matrix_pre = cov(normalised_pre');
[V, E] = eig(cov_matrix_pre);
%gives eigenvectors as columns in V and diagonal matrix E of eigen values 
[eigen_values, columns] = sort(diag(E),'descend');
%sorts eigen values in descending and gives columns indices
E_sorted = E(columns, columns);
%sorts E in descending order (PC1 to PC85)
V_sorted = V(:, columns);
%sorts E in descending order of columns (PC1 to PC85)

cov_test_score = V_sorted * normalised_pre;
cov_test_scores_post = V_sorted * normalised_post;

for i = 1:numel(eigen_values)
    [~, eigenval_col(i)] = find(E == eigen_values(i));
    eigenvec_col(i) = V(:, eigenval_col(i));
end
%gives the transforming vector for PC1

sum_coeff_1_sq = sum(coeff_1.^2);
recreate_score =  normalised_pre * u;
recreate_score_no_mean = coeff_1 * pre_minus_mean;
difference = recreate_score_no_mean - score_1;

scores_post = coeff_1 * normalised_post;
%coeff_1 each column is a PC, each row is the weighting of each neuron, the
%coefficient for that neuron in a vector of all neurons that correspond to
%the eigenvector of the new orthonormal basis

%plots coefficient for all neurons for 1st 2 PCs, 
c = turbo(nframe);

%prestim
figure;hold on;
title('PCA scores, prestim');
color_bar = colorbar;
bar_feature_1 = color_bar.Limits;
color_bar.Limits = [0 nframe];
for i = 1:nframe;
    %biplot(coeff(:,1:2),'scores',score(:,1:2),)
    plot(score_1(i,1), score_1(i,2), '.', 'color', c(i,:))
end
for i = 1:npair;
    plot([0, coeff_1(i,1)],[0, coeff_1(i,2)], 'color', 'b')
    text(coeff_1(i,1), coeff_1(i,2), num2str(i))
end

figure; hold on;
for cc = 1:nframe;
    plot(cc, cc, '.', 'color', c(cc,:));
end
%poststim
figure;hold on;
title('PCA scores, poststim');
for i = 1:nframe;
    %biplot(coeff(:,1:2),'scores',score(:,1:2),)
    plot(score_2(i,1), score_2(i,2), '.', 'color', c(i,:))
end
for i = 1:npair;
    plot([0, coeff_2(i,1)],[0, coeff_2(i,2)], 'color', 'b')
    text(coeff_2(i,1), coeff_2(i,2), num2str(i))
end

%plot with no movement





%for the biplot one
for i = 1:numel(pre.F(:,1))    
    a{i} = num2str(i)
end
c = turbo(nframe);
figure;hold on;
for i = 1:nframe;
    biplot(coeff_1(:,1:2),'scores',score_1(:,1:2));
end


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
