% function align_images
%% load the figs and plotting
[figname1,path1a]=uigetfile('.fig','find the fig');

fig1=openfig(fullfile(path1a,figname1));
fig1img =getimage(fig1);
kC1= fig1.Children.Children;
cellnum1=(numel(kC1)-1)/2;

SIZ = size(fig1img);
F1Contour =zeros(prod(SIZ),cellnum1);

for j=1:cellnum1
    nth = cellnum1+1-j;
    F1cell(j).Position = [kC1(nth).Position(1)+v+5,kC1(nth).Position(2)+u+5];
    F1cell(j).String = kC1(nth).String;
    F1Contour(:,j) = reshape(circshift(kC1(nth+cellnum1).ZData,[u,v]),prod(SIZ),1);
end

pos=strfind(figname1,'_');
name=[figname1(1:pos(4)-1) '&' figname2(1:pos(4)-1) ];

%% Select good cells and plot
prompt = 'red cell';
title = 'input the green cell# matching the red';
pair1=[];
pair2=[];

h1=figure;
imshow(fig1img)
set(gcf,'Position',[0 0 2000 1200]);
hold on;

%2nd image, green
for j=1:cellnum1
        text(F1cell(j).Position(1),F1cell(j).Position(2),F1cell(j).String,'Color',[ 1 0 0],'FontSize',8)
        contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,j),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[1 0 0])
        drawnow;
        default ='';
        
        gr=inputdlg({[prompt num2str(j)]}, title ,[2 80],{default});
        if ~isempty(gr{1})
            i2=str2num(gr{1});
            temp= F1Contour(:,j).*F2Contour(:,i2);
            if sum(temp(:)) == 0
                disp('wrong pairs? Redo matching'); j = j-1;
            else
                pair1=[pair1;j];
                text(F1cell(j).Position(1),F1cell(j).Position(2),F1cell(j).String,'Color',[ 1 1 1],'FontSize',8)
                contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,j),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[1 1 1])
               
            end
        end
    end
    j = j+1;
end

%%
fig4=figure;
imshowpair(fig1imgnew,fig2img)
set(gcf,'Position',[0 0 2000 1200]);
hold on;
%first image,red
for m=1:numel(pair1)
    text(F1cell(pair1(m)).Position(1),F1cell(pair1(m)).Position(2),F1cell(pair1(m)).String,'Color',[ 1 0 0],'FontSize',8)
    contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,pair1(m)),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[1 0 0])
end
%% 
saveas(fig4,[name '_selected'],'fig');
saveas(fig4,[name '_selected'],'png');
save([name '_selected.mat' ],'pair1','pair2');

close all