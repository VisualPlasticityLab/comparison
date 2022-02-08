% Author: Aashish Khimasia
% Date: July 2021
%% Getting Fall.mat files and storing in tmp
for ii=1:3 
    % ii is number of the plane (1 to 3)
    tmp{ii}=load(fullfile(['plane' num2str(ii-1)],'Fall.mat'));   
end
%% Numbers of cells for each plane to be plotted
cell_num = [[4 104 140], 
    [19 38 71], 
    [17 40 4]];
%% Plot figures 
for b = 1:3;
    figure;hold on;
    for a = 1:3
        subplot(3,1,a);
        cell = cell_num(b, a);
        plot(tmp{b}.F(cell,:), color = 'b'); hold on;
        plot(tmp{b}.Fneu(cell,:), color = 'r'); hold on;
        c = turbo(4);
        for i = 7:10
            coeff = i/10;
            plot(tmp{b}.F(cell,:) - (coeff * tmp{b}.Fneu(cell, :)), 'color', c(i-6, :))
            hold on;
        end
        yline(0);
        legend({'F', 'Fneu', 'F-0.7*Fneu', 'F-0.8*Fneu','F-0.9*Fneu','F-Fneu'}, 'Orientation', 'horizontal');
        figaxes(a) = gca;
        figaxes(a).XLim = [0 2400];
        figaxes(a).YLim = [min(tmp{b}.F(cell,:))-1000 max(tmp{b}.F(cell,:))+2000];
        figaxes(a).XTick = [0 400 800 1200 1600 2000 2400];
        figaxes(a).XLabel.String = 'Frame number';
        figaxes(a).YLabel.String = 'Activity';
        figaxes(a).Title.String = ['Cell ' num2str(cell) ' Plane ' num2str(b-1)];
    end
end