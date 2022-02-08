% plot skewness for redcell and nonredcells seperately
function plot_skewness_rc(pre, post, redcell)
% redcell_index = find(post.redcell2(:,4)==1);
% non_redcell_index = find(post.redcell2(:,4)==0);

for m = 1:5
    redcell_index{m} = find(redcell{m}(:,1)==1);
    non_redcell_index{m} = find(redcell{m}(:,1)==0);
    for i = 1:numel(redcell_index{m});
        red_skew_pre{m}(i) = pre{m}.stat{redcell_index{m}(i)}.skew;
        red_skew_post{m}(i) = post{m}.stat{redcell_index{m}(i)}.skew;
    end

    for j = 1:numel(non_redcell_index{m});
        non_red_skew_pre{m}(j) = pre{m}.stat{non_redcell_index{m}(j)}.skew;
        non_red_skew_post{m}(j) = post{m}.stat{non_redcell_index{m}(j)}.skew;
    end
    
end
    all_red_skew_pre = [red_skew_pre{1},red_skew_pre{2},red_skew_pre{3},red_skew_pre{4},red_skew_pre{5}];
    all_red_skew_post = [red_skew_post{1},red_skew_post{2},red_skew_post{3},red_skew_post{4},red_skew_post{5}];
    all_non_red_skew_pre = [non_red_skew_pre{1},non_red_skew_pre{2},non_red_skew_pre{3},non_red_skew_pre{4},non_red_skew_pre{5}];
    all_non_red_skew_post = [non_red_skew_post{1},non_red_skew_post{2},non_red_skew_post{3},non_red_skew_post{4},non_red_skew_post{5}];

%% Plots coloured box plots with error bars
figure;hold on;
title('Skewness')
positions_1 = [1;2];
positions_2 = [3.5;4.5];
mean_data_red = [mean(all_red_skew_pre); mean(all_red_skew_post)];
mean_data_non_red = [mean(all_non_red_skew_pre); mean(all_non_red_skew_post)];
std_data_red = [std(all_red_skew_pre); std(all_red_skew_pre)];
std_data_non_red = [std(all_non_red_skew_post); std(all_non_red_skew_post)];
plot([1.1;2.1],[all_red_skew_pre', all_red_skew_post'],'x-', 'Color',[.5 .5 .5 .3], 'LineWidth',.5)
errorbar(positions_1, mean_data_red, std_data_red, 'O', 'color', 'r', 'LineWidth', 2);
plot([3.6;4.6],[all_non_red_skew_pre', all_non_red_skew_post'],'x-', 'Color',[.5 .5 .5 .3], 'LineWidth',.5)
errorbar(positions_2, mean_data_non_red, std_data_non_red, 'O', 'color', 'k', 'LineWidth', 2);
figaxes(1) = gca;
figaxes(1).Title.String = ['Skewness, all mice'];
figaxes(1).XLabel.String = 'cells';
figaxes(1).YLabel.String = 'skewness';
figaxes(1).XTick = [1 2 3.5 4.5];
figaxes(1).XLim = [0.5 5];
figaxes(1).XTickLabel = ({'pre-red', 'post-red', 'pre-nonred', 'post-nonred'});
%% Plots separate box plots
% figure; hold on;
% boxplot([red_skew_pre', red_skew_post'])
% plot([1;2],[red_skew_pre', red_skew_post'],'x-', 'Color',[.5 .5 .5 .3], 'LineWidth',.5)
% title([pre.ops.mouse_name ' redcell skewness, pre and post'])
% xlim([.5 2.5]);
% xticks(1:2);
% xticklabels({'pre-stimulation-skew','post-stimulation-skew'})
% 
% figure; hold on;
% boxplot([non_red_skew_pre', non_red_skew_post'])
% plot([1;2],[non_red_skew_pre', non_red_skew_post'],'x-', 'Color',[.5 .5 .5 .3], 'LineWidth',.5)
% title([pre.ops.mouse_name ' non redcell skewness, pre and post'])
% xlim([.5 2.5]);
% xticks(1:2);
% xticklabels({'pre-stimulation-skew','post-stimulation-skew'})

%% Plots error bars
% figure;hold on;
% title('Skewness')
% subplot(1,2,1);
% mean_red = [mean(red_skew_pre); mean(red_skew_post)];
% std_red = [std(red_skew_pre); std(red_skew_post)];
% errorbar(mean_red, std_red, 'O');
% figaxes(2) = gca;
% figaxes(2).Title.String = 'Skewness of red cells (pre post)';
% figaxes(2).XLabel.String = 'red cells';
% figaxes(2).YLabel.String = 'skewness';
% figaxes(2).XTick = [1 2];
% figaxes(2).XLim = [0.5 2.5];
% figaxes(2).XTickLabel = ({'pre', 'post'});
% 
% % errorbar(mean(red_skew_pre), std(red_skew_pre), 'x')
% % errorbar(mean(red_skew_post), std(red_skew_post), '.')
% 
% subplot(1,2,2);
% mean_non_red = [mean(non_red_skew_pre); mean(non_red_skew_post)];
% std_non_red = [std(non_red_skew_pre); std(non_red_skew_post)];
% errorbar(mean_non_red, std_non_red, 'O');
% figaxes(3) = gca;
% figaxes(3).Title.String = 'Skewness of non red cells (pre post)'
% figaxes(3).XLabel.String = 'non red cells';
% figaxes(3).YLabel.String = 'skewness';
% figaxes(3).XLim = [0.5 2.5];
% figaxes(3).XTick = [1 2];
% figaxes(3).XTickLabel = ({'pre', 'post'});
% % errorbar(mean(non_red_skew_pre), std(non_red_skew_pre), 'x')
% % errorbar(mean(non_red_skew_post), std(non_red_skew_post), '.')