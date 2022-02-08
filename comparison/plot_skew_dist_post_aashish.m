%Plotting distance from electrode against skewness

%Issues = what plane to we use. Do we use pre and post? Would they have rel
%location data? We need to get skew and rel loc from the same structure
%pre and post have med, but not the same. Are we doing change in skewness?
%get from Fall.mat
%change of skewness
%post position, post - pre change of skewness

plane2ak = cell(5, 1);

plane2ak{1} = './GAD1/190214/002/suite2p/plane2/Fall.mat';
%plane2ak{2} = './GADA/190128/002/suite2p/plane2/Fall.mat';
plane2ak{2} = './GADB/190128/004/suite2p/plane2/Fall.mat';
plane2ak{3} = './GADD/190128/002/suite2p/plane2/Fall.mat';
plane2ak{4} = './GADC/190130/002/suite2p/plane2/tiffs/Fall.mat';

section_name = strsplit(plane2ak{1}, '/');
mouse_name = section_name(2);
epos = load('Z:\Teams-MJ-data\epos.mat');

for i = 1:numel(plane2ak)
    curr_plane_data = load(plane2ak{i});
    iscell = cat(1,curr_plane_data.iscell); 
    cell_indices{i} = find(iscell(:,1) == 1);
    for j = 1:numel(cell_indices{i})
        temp_loc{j} = curr_plane_data.stat{j}.med;
        temp_loc_test(j,1) = curr_plane_data.stat{j}.med(1,1);
        temp_loc_test(j,2) = curr_plane_data.stat{j}.med(1,2);

%         tmplocs = cat(1,curr_plane_data.stat{j}.med);
        temp_skew(j) = curr_plane_data.stat{j}.skew;
%         skew = cat(1,curr_plane_data.stat.skew);
%     cell_locs = temp_loc
    cell_skew{i} = cell2mat(cat(1, temp_skew))
%         cell_locs = tmplocs(cell_indices{i}, :);
    rel_locs = cell_locs - epos.epos(i, 1:2);
    cell_dist{i} = sqrt(rel_locs(:,1).^2 + rel_locs(:,2).^2);
    cell_skew{i} = skew(cell_indices{i}, :);
end

figure; hold on;
title('Cell skewness vs distance');
for m = 1:numel(plane2ak)
    subplot(2,3,m); hold on;
    for j = 1:numel(cell_indices{m})
        plot(cell_dist{m}(j), cell_skew{m}(j), 'x', 'color', 'b')
    figaxes(m) = gca;
    figaxes(m).XLabel.String = 'distance from electrode um';
    figaxes(m).YLabel.String = 'skewness';
    figaxes(m).Title.String = [mouse_name 'dist vs skew'];
    end
end


%FOR MARIAS F MAT FILES
% for i = 1:numel(plane2ak)
%     curr_plane_data = load(plane2ak{i});
%     iscell = cat(1,curr_plane_data.stat.iscell); 
%     cell_indices{i} = find(iscell(:,1) == 1);
%     tmplocs = cat(1,curr_plane_data.stat.med);
%     skew = cat(1,curr_plane_data.stat.skew);
%     cell_locs = tmplocs(cell_indices{i}, :);
%     rel_locs = cell_locs - epos.epos(i, 1:2);
%     cell_dist{i} = sqrt(rel_locs(:,1).^2 + rel_locs(:,2).^2);
%     cell_skew{i} = skew(cell_indices{i}, :);
% end
