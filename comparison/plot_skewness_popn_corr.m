%% Skewness vs Popn correlation
function plot_skewness_popn_corr(pre, post, movie_frames_pre, movie_frames_post, redcell, npair)

for m = 1:5
    normalised_pre{m} = zscore(pre{m}.F(:,movie_frames_pre{m}), 0, 2);
    normalised_post{m} = zscore(post{m}.F(:,movie_frames_post{m}), 0, 2);

    %change below for pre or post
    popn_avg{m} = mean(normalised_pre{m});
    % popn_avg = mean(normalised_post);
    for i = 1:npair{m}
        R{m} = corrcoef(popn_avg{m}, normalised_pre{m}(i,:));
    %   R = corrcoef(popn_avg, normalised_post(i,:));
        popn_corr_coeff{m}(i) = R{m}(1,2);
        skewness_change{m}(i) = post{m}.stat{i}.skew - pre{m}.stat{i}.skew;
    end
    
    redcell_index{m} = find(redcell{m}(:,1) == 1);
    non_redcell_index{m} = find(redcell{m}(:,1) == 0);
    
    skewness_change_rc{m} = skewness_change{m}(redcell_index{m});
    skewness_change_nonrc{m} = skewness_change{m}(non_redcell_index{m});

    popn_corr_rc{m} = popn_corr_coeff{m}(redcell_index{m});
    popn_corr_nonrc{m} = popn_corr_coeff{m}(non_redcell_index{m});

end

    all_skew_rc = [skewness_change_rc{1},skewness_change_rc{2},skewness_change_rc{3},skewness_change_rc{4},skewness_change_rc{5}];
    all_skew_nonrc = [skewness_change_nonrc{1},skewness_change_nonrc{2},skewness_change_nonrc{3},skewness_change_nonrc{4},skewness_change_nonrc{5}];
    all_popncorr_rc = [popn_corr_rc{1},popn_corr_rc{2},popn_corr_rc{3},popn_corr_rc{4},popn_corr_rc{5}];
    all_popncorr_nonrc = [popn_corr_nonrc{1},popn_corr_nonrc{2},popn_corr_nonrc{3},popn_corr_nonrc{4},popn_corr_nonrc{5}];

    figure; hold on;
    title('Skewness change vs popn correlation, all mice')
    %   title([mouse_name ' skewness change vs popn correlation (post)'])
    plot(all_skew_rc, all_popncorr_rc, '.', 'MarkerSize', 7, 'color', 'r');
    plot(all_skew_nonrc, all_popncorr_nonrc, 'x', 'MarkerSize', 5, 'color', 'k');      
    legend('inhibitory', 'excitatory')
    xlabel('skewness change (post - pre)')
    ylabel('population correlation')
    end
    
%     all_skewness_change = [skewness_change{1},skewness_change{2},skewness_change{3},skewness_change{4},skewness_change{5}];
%     skewness_change_red = skewness_change(redcell_index);
%     skewness_change_nonred = skewness_change(non_redcell_index);
%     popn_corr_red = popn_corr_coeff(redcell_index);
%     popn_corr_nonred = popn_corr_coeff(non_redcell_index);
% 
%     % Linear regression
%     [p, S] = polyfit(skewness_change, popn_corr_coeff, 1);
%     line_of_best_fit = polyval(p, skewness_change);
%     plot(skewness_change, line_of_best_fit, '-x');
% 
%     [p_red, S_red] = polyfit(skewness_change_red, popn_corr_red, 1);
%     [p_nonred, S_nonred] = polyfit(skewness_change_nonred, popn_corr_nonred, 1);
%     line_of_best_fit_red = polyval(p_red, skewness_change_red);
%     line_of_best_fit_nonred = polyval(p_nonred, skewness_change_nonred);
%     plot(skewness_change_red, line_of_best_fit_red, '-x', 'color', 'r');
%     plot(skewness_change_nonred, line_of_best_fit_nonred, '-x', 'color', 'k');