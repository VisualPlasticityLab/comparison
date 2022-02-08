% Author = Aashish Khimasia
% Date = July 2021
%% create list of duplicates from green cells and indices of corresponding red pairs
green_duplicates = [];
red_indices = [];
sorted_pair2 = sort(pair2);
for i = 1:(numel(pair2) - 1)
    if sorted_pair2(i) == sorted_pair2(i+1);
        green_duplicates(end+1) = sorted_pair2(i);
    end
end
green_duplicates
for i = 1:numel(green_duplicates)
    red_indices{end+1} = find(pair2==green_duplicates(i));
end
red_indices


%% visualise the duplicates on the figure image
for j = 1:numel(green_duplicates)
    figtemp=figure;
    fig1imgnew = circshift(fig1img,[u,v]);
    imshowpair(fig1imgnew,fig2img)
    set(gcf,'Position',[0 0 1200 1000]);
    hold on;
    %first image,red 
    for m=[red_indices{j}]
        text(F1cell(pair1(m)).Position(1),F1cell(pair1(m)).Position(2),F1cell(pair1(m)).String,'Color',[ 1 0 0],'FontSize',8)
        contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,pair1(m)),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[1 0 0])

        text(F2cell(pair2(m)).Position(1),F2cell(pair2(m)).Position(2),F2cell(pair2(m)).String,'Color',[ 0 1 0],'FontSize',8)
        contour(1:SIZ(2),1:SIZ(1),reshape(F2Contour(:,pair2(m)),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[0 1 0])

    end
    
end

%% find the numbers of the remaining cells with low overlap

low_overlap = [];
for i = 1:numel(low_overlap_temp)
    for j = 1:numel(pair1)
        if pair1(j) == low_overlap_temp(i)
            low_overlap(end+1) = low_overlap_temp(i);
        end
    end
end

%% print image with only poor overlap cells
figtemp1=figure;
title('poorly matched cells')
fig1imgnew = circshift(fig1img,[u,v]);
imshowpair(fig1imgnew,fig2img)
set(gcf,'Position',[0 0 1200 1000]);
hold on;
%first image,red
for m=1:numel(low_overlap)
    n = find(pair1==low_overlap(m));
    text(F1cell(pair1(n)).Position(1),F1cell(pair1(n)).Position(2),F1cell(pair1(n)).String,'Color',[ 1 0 0],'FontSize',8)
    contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,pair1(n)),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[1 0 0])
    
    text(F2cell(pair2(n)).Position(1),F2cell(pair2(n)).Position(2),F2cell(pair2(n)).String,'Color',[ 0 1 0],'FontSize',8)
    contour(1:SIZ(2),1:SIZ(1),reshape(F2Contour(:,pair2(n)),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[0 1 0])
    
end