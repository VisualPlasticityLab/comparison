%% load all pairs data
se = load('.\F_GAD1_allpairs_JA.mat');
%% create pre and post
% function [pre, post, neur_coeff, redcell, npair, nframe] = create_pre_post(se)
neur_coeff = 0.7;
pre.F =se.F_pre - neur_coeff*se.Fneu_pre;
post.F =se.F_post - neur_coeff*se.Fneu_post;
pre.stat = se.stat_pre; 
post.stat = se.stat_post; 
pre.ops = se.ops_pre{1};
post.ops = se.ops_post{1};
redcell = se.redcell;
pre.velocity = se.velocity_pre;
post.velocity = se.velocity_post;

npair = size(pre.F,1);
nframe = size(pre.F,2); 
% end
%% set movie frames
% function [movie_frames_pre, movie_frames_post, movie_sect_length] = set_movie_frames(pre, post, nframe) 
movie_frames_pre = [];
movie_frames_post = [];

all_start_frame_pre = ceil((pre.ops.start_frame-1)/3+1);
all_start_frame_post = ceil((post.ops.start_frame-1)/3+1);

if all_start_frame_pre >= all_start_frame_post
    start_frame_later = all_start_frame_pre;
else
    start_frame_later = all_start_frame_post;
end

movie_sect_length = nframe - start_frame_later + 1;

start_frame_pre = all_start_frame_pre;
start_frame_post = all_start_frame_post;
end_frame_pre = all_start_frame_pre + movie_sect_length - 1;
end_frame_post = all_start_frame_post + movie_sect_length - 1;

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
% end

%% plot skewness
% function plot_skewness(pre, post)
for i = 1:numel(pre.stat);
    pre_skew(i) = pre.stat{i}.skew;
end
for j = 1:numel(post.stat);
    post_skew(j) = post.stat{j}.skew;
end

figure; hold on;
boxplot([pre_skew', post_skew']);
plot([1;2], [pre_skew', post_skew'],'x-', 'Color',[.5 .5 .5 .3], 'LineWidth',.5);
title([pre.ops.mouse_name ' skewness, pre and post'])
xlim([.5 2.5]);
xticks(1:2);
xticklabels({'pre-stimulation-skew','post-stimulation-skew'})
% end
%% Skewness vs Popn correlation
% function = plot_skew_change_popn_corr(pre, post, movie_frames_pre,
% movie_frames_post, redcell)
normalised_pre = zscore(pre.F(:,movie_frames_pre), 0, 2);
normalised_post = zscore(post.F(:,movie_frames_post), 0, 2);

%change below for pre or post
popn_avg = mean(normalised_pre);
% popn_avg = mean(normalised_post);
for i = 1:npair
    R = corrcoef(popn_avg, normalised_pre(i,:));
%     R = corrcoef(popn_avg, normalised_post(i,:));
    popn_corr_coeff(i) = R(1,2);
    skewness_change(i) = post.stat{i}.skew - pre.stat{i}.skew;
end

redcell_index = find(redcell(:,1) == 1);
non_redcell_index = find(redcell(:,1) == 0);
figure; hold on;
title([pre.ops.mouse_name ' skewness change vs popn correlation (pre)'])
% title([mouse_name ' skewness change vs popn correlation (post)'])
c = turbo(npair);
for j = 1:numel(redcell_index)
    k = redcell_index(j);
    plot(skewness_change(k), popn_corr_coeff(k), '.', 'MarkerSize', 7, 'color', 'r');
end
for h = 1:numel(non_redcell_index)
    m = non_redcell_index(h);
    plot(skewness_change(m), popn_corr_coeff(m), 'x', 'MarkerSize', 5, 'color', 'k');
end
legend('inhibitory', 'excitatory')
xlabel('skewness change (post - pre)')
ylabel('population correlation')
% color_bar = colorbar;
% color_bar.Ticks = [0 1];
% color_bar.TickLabels = {'1', num2str(npair)};
% colormap turbo;

% Linear regression
[p, S] = polyfit(skewness_change, popn_corr_coeff, 1);
% line(skewness_change, p(1)*skewness_change + p(2));
line_of_best_fit = polyval(p, skewness_change)
plot(skewness_change, line_of_best_fit, '-x')
% end
%% pca analysis
function [score] = plot_pca(pre, post, movie_frames_pre, movie_frames_post, movie_sect_length)

%normalising the data
normalised_pre = zscore(pre.F(:,movie_frames_pre), 0, 2);
normalised_post = zscore(post.F(:,movie_frames_post), 0, 2);
normalised = cat(2,normalised_pre,normalised_post);

%run pca analysis
[coeff, score, latent, ~, explained] = pca(normalised');

%eigenvectors, scores(remapping of data), eigenvalues, ~, percentage of
%variance explained by the PC

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
figaxes(1).Title.String = [pre.ops.mouse_name ' pre-stim'] % thr_status];
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
figaxes(2).Title.String = [post.ops.mouse_name ' post-stim'] % thr_status];
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

end
%% calculate pca correlation between pre and post
function plot_pca_corr(pre, post, movie_sect_length, score)

for i = 1:size(pre.F,1)
    R = corrcoef(score(1:movie_sect_length,i), score((movie_sect_length + 1):movie_sect_length * 2,i));
    pca_corr_coeff(i) = R(1,2);
end

figure;hold on;
for i = 1:size(pre.F,1)
    plot(i, pca_corr_coeff(i), 'x', 'color', 'r')
end
title([pre.ops.mouse_name ': corr pre vs post for each PC']);
xlabel('PC');
ylabel('correlation_coefficient');
end