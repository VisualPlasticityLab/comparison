%get skewness information
%plot this in a fig

for i = 1:numel(pre.stat);
    pre_skew(i) = pre.stat{i}.skew;
end
for j = 1:numel(post.stat);
    post_skew(j) = post.stat{j}.skew;
end
file_name = strsplit(se2{1}.ops.data_path, '/');
mouse_name = file_name(3);
figure; hold on;
boxplot([pre_skew', post_skew']);
plot([1;2], [pre_skew', post_skew'],'x-', 'Color',[.5 .5 .5 .3], 'LineWidth',.5);
title([mouse_name 'skewness, pre and post'])
xlim([.5 2.5]);
xticks(1:2);
xticklabels({'pre-stimulation-skew','post-stimulation-skew'})

