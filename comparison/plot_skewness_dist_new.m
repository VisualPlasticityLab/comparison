function [rellocs] = plot_skewness_dist_prepost(m, pre, post, npair, epos, redcell)
for i = 1:npair

    skewness_change(i) = post.stat{i}.skew - pre.stat{i}.skew;
    tmplocs(i,1) = pre.stat{i}.med(1);
    tmplocs(i,2) = pre.stat{i}.med(2);
    tmplocplanes(i) = pre.stat{i}.iplane;

end
tmplocs = double(tmplocs);
locs = tmplocs.* repmat(epos.micronppixel/epos.mag(m),size(tmplocs(:,1))); 
rellocs = (locs - epos.epos(m,1:2)) .* repmat(epos.micronppixel/epos.mag(m),size(locs,1),1); 
locs(:,3) = tmplocplanes * epos.micronpplane;
rellocs(:,3)=(tmplocplanes - epos.epos(m,3))*epos.micronpplane;

dist = sqrt(rellocs(:,1).^2 + rellocs(:,2).^2+rellocs(:,3).^2);

redcell_index = find(redcell(:,1) == 1);
non_redcell_index = find(redcell(:,1) == 0);

figure; hold on;
title([pre.ops.mouse_name ' skewness change vs distance from electrode'])

for j = 1:numel(redcell_index)
    k = redcell_index(j);
    plot(skewness_change(k), dist(k), '.', 'MarkerSize', 7, 'color', 'r');
end
for h = 1:numel(non_redcell_index)
    n = non_redcell_index(h);
    plot(skewness_change(n), dist(n), 'x', 'MarkerSize', 5, 'color', 'k');
end
legend('inhibitory', 'excitatory')
xlabel('skewness change (post - pre)')
ylabel('distance')


skewness_change_red = skewness_change(redcell_index);
skewness_change_nonred = skewness_change(non_redcell_index);
dist_red = dist(redcell_index);
dist_nonred = dist(non_redcell_index);
% Linear regression
[p, S] = polyfit(skewness_change, dist, 1);
line_of_best_fit = polyval(p, skewness_change);
plot(skewness_change, line_of_best_fit, '-x');

[p_red, S_red] = polyfit(skewness_change_red, dist_red, 1);
[p_nonred, S_nonred] = polyfit(skewness_change_nonred, dist_nonred, 1);
line_of_best_fit_red = polyval(p_red, skewness_change_red);
line_of_best_fit_nonred = polyval(p_nonred, skewness_change_nonred);
plot(skewness_change_red, line_of_best_fit_red, '-x', 'color', 'r');
plot(skewness_change_nonred, line_of_best_fit_nonred, '-x', 'color', 'k');
