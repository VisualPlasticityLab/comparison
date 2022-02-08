function [rellocs] = plot_skewness_dist_prepost(pre, post, npair, epos, redcell)
for m = 1:5
    for i = 1:npair{m}

        skewness_change{m}(i) = post{m}.stat{i}.skew - pre{m}.stat{i}.skew;
        tmplocs{m}(i,1) = pre{m}.stat{i}.med(1);
        tmplocs{m}(i,2) = pre{m}.stat{i}.med(2);
        tmplocplanes{m}(i) = pre{m}.stat{i}.iplane;

    end
    tmplocs{m} = double(tmplocs{m});
    locs{m} = tmplocs{m}.* repmat(epos.micronppixel/epos.mag(m),size(tmplocs{m}(:,1))); 
    rellocs{m} = (locs{m} - epos.epos(m,1:2)) .* repmat(epos.micronppixel/epos.mag(m),size(locs{m},1),1); 
    locs{m}(:,3) = tmplocplanes{m} * epos.micronpplane;
    rellocs{m}(:,3)=(tmplocplanes{m} - epos.epos(m,3))*epos.micronpplane;

    dist{m} = sqrt(rellocs{m}(:,1).^2 + rellocs{m}(:,2).^2+rellocs{m}(:,3).^2);

    redcell_index{m} = find(redcell{m}(:,1) == 1);
    non_redcell_index{m} = find(redcell{m}(:,1) == 0);
    
    skewness_change_rc{m} = skewness_change{m}(redcell_index{m});
    skewness_change_nonrc{m} = skewness_change{m}(non_redcell_index{m});
    dist_rc{m} = dist{m}(redcell_index{m});
    dist_nonrc{m} = dist{m}(non_redcell_index{m});

end
    
    all_skewness_change=[skewness_change{1},skewness_change{2},skewness_change{3},skewness_change{4},skewness_change{5}];
    all_dist = [dist{1};dist{2};dist{3};dist{4};dist{5}];
    all_skew_rc = [skewness_change_rc{1},skewness_change_rc{2},skewness_change_rc{3},skewness_change_rc{4},skewness_change_rc{5}];
    all_skew_nonrc = [skewness_change_nonrc{1},skewness_change_nonrc{2},skewness_change_nonrc{3},skewness_change_nonrc{4},skewness_change_nonrc{5}];
    all_dist_rc = [dist_rc{1};dist_rc{2};dist_rc{3};dist_rc{4};dist_rc{5}];
    all_dist_nonrc = [dist_nonrc{1};dist_nonrc{2};dist_nonrc{3};dist_nonrc{4};dist_nonrc{5}];

    figure; hold on;
    plot(all_dist_rc, all_skew_rc, '.', 'MarkerSize', 7, 'color', 'r');
    plot(all_dist_nonrc, all_skew_nonrc,  'x', 'MarkerSize', 7, 'color', 'k');

    legend('inhibitory', 'excitatory')
    ylabel('skewness change (post - pre)')
    xlabel('distance (um)')
    title('Skewness change vs distance from electrode, all mice')

% figure; hold on;
% title('Skewness change vs distance from electrode, all mice')
% for m = 1:5
%     for j = 1:numel(redcell_index{m})
%         k = redcell_index{m}(j);
%         plot(dist{m}(k), skewness_change{m}(k), '.', 'MarkerSize', 7, 'color', 'r');
%     end
%     for h = 1:numel(non_redcell_index{m})
%         n = non_redcell_index{m}(h);
%         plot(dist{m}(n), skewness_change{m}(n), 'x', 'MarkerSize', 5, 'color', 'k');
%     end
%     legend('inhibitory', 'excitatory')
%     ylabel('skewness change (post - pre)')
%     xlabel('distance (um)')
% end
%% 
% db = 50;
% distlist = [0:db:500, 600];
% distcent = distlist(2)/2 + distlist(1:end-1);
% skewness_change_red = skewness_change(redcell_index);
% skewness_change_nonred = skewness_change(non_redcell_index);
% dist_red = dist(redcell_index);
% dist_nonred = dist(non_redcell_index);
% % Linear regression
% [p, S] = polyfit(skewness_change, dist, 1);
% line_of_best_fit = polyval(p, skewness_change);
% plot(skewness_change, line_of_best_fit, '-x');
% 
% [p_red, S_red] = polyfit(skewness_change_red, dist_red, 1);
% [p_nonred, S_nonred] = polyfit(skewness_change_nonred, dist_nonred, 1);
% line_of_best_fit_red = polyval(p_red, skewness_change_red);
% line_of_best_fit_nonred = polyval(p_nonred, skewness_change_nonred);
% plot(skewness_change_red, line_of_best_fit_red, '-x', 'color', 'r');
% plot(skewness_change_nonred, line_of_best_fit_nonred, '-x', 'color', 'k');

%%
% figure(103)
% colorlist = colorbrewerRGB(nampl, 'qualitative');
% for k = 1:nampl
% %     subplot(2,2,1), hold on
% %      errorbar(distlist(1:end-1) + db/2, meand(k,:), stdd(k,:),'o-','linewidth',2,'color',colorlist(k,:))
%     subplot(2,1,1), hold on
%     errorbar(distcent(binlist), meancom(k,:), stdcom(k,:),'o-','linewidth',2,'color',colorlist(k,:))
% end
% xlabel('Cell Distance')
% ylabel('Average Evoked Response')
% title('Cell activation as a function of Stimulation Amplitude and Cell Location COM')
% legend(num2str(amplid))
% cc = turbo(nbin);
% for g = 1:nbin
% %     subplot(2,2,3), hold on
% %     errorbar(amplid, meand(:,g), stdd(:,g),'o-','linewidth',2,'color',cc(g,:))
%     subplot(2,1,2), hold on
%     errorbar(amplid, meancom(:,g), stdcom(:,g),'o-','linewidth',2,'color',cc(g,:))
% end
% xlabel('Stimulation Amplitude')
% ylabel('Average Evoked Response')
% legend(num2str(distlist(:)+db/2))
% title('Cell activation as a function of Stimulation Amplitude and Cell Location COM')
% set(gcf,'Position',[10 10 500 700])
% saveas(gcf,['Distr_of_Resps_Ampl_Dist_COM_',TraceType,'.pdf'])
