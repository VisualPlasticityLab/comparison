% function align_two_days
% this function compare two images and find physically matched cell, check
% cell responses with the function analysis_twoday_data, see more in Demo
% input: to be selected figures
% output: automatically saved matfile 
%% (redo step 1)load the files and figs 
[figname1,path1a]=uigetfile('.fig','find the figure to be compared');
[figname2,path2a]=uigetfile('.fig','find the 2nd baseline fig');
pos=strfind(strtok(figname1,'&'),'_');
pos_plane1 = strfind(figname1,'plane');
pos_plane2 = strfind(figname2,'plane');
nplane = figname1(pos_plane1+5);
name=[figname1(1:pos_plane1-2) '&&' figname2(1:pos_plane2-2) '_plane' nplane];

fig1=openfig(fullfile(path1a,figname1));
fig1img =getimage(fig1);
% fig1img = squeeze(fig1img(:,:,1)); %red channel
if ndims(fig1img)>2
    fig1img = fig1img(:,:,2); %green channel 
end
    kC1= fig1.Children.Children;
cellnum1=(numel(kC1)-1)/2;

fig2=openfig(fullfile(path2a,figname2));
fig2img= getimage(fig2);
% fig2img = squeeze(fig2img(:,:,1)); %red channel
if ndims(fig2img)>2
    fig2img = fig2img(:,:,2); %green channel
end
kC2= fig2.Children.Children;
cellnum2=(numel(kC2)-1)/2;
%% align two images and plot figure 1 in red & figure2 in green
fig3=openfig(fullfile(path2a,figname2));
while menu('Select areas for image alignment?','Yes','Done')==1
    rect=round(getrect(fig3));
    img1=fig1img(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3));
    img2=fig2img(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3));
    %     max1 = max(img1(:));
    %     min1 = min(img1(:));
    %     max2 = max(img2(:));
    %     min2 = min(img2(:));
    %     img1=(img1-min1)/(max1-min1)*(max2-min2)+min2;
    %     [optimizer, metric] = imregconfig('multimodal');
    %     moving= img1;
    %     fixed = img2;
    %     [moving_reg,R_reg] =  imregister(moving, fixed, 'rigid', optimizer, metric);
    %     tform = imregtform(moving, fixed, 'rigid', optimizer, metric);
    %     v=round(tform.T(3,1));
    %     u=round(tform.T(3,2));
    [u,v] = fftalign(img1,img2);
    fprintf('registration shift(%d,%d)\n',u,v);
    fig1imgnew = circshift(fig1img,[u,v]);
    figure(fig3);
    imshowpair(fig1imgnew,fig2img)
    title('fig1 in red & fig2 in green')
%     set(fig3,'Position',[0 0 2000 1200]);
end
%% (redo step 2) read the cell position and text info
SIZ = size(fig1img);

F1Contour =zeros(prod(SIZ),cellnum1);
F2Contour =zeros(prod(SIZ),cellnum2);
Paninski = strcmp(class(kC1(1)),class(kC1(2)));
if Paninski
    %1st image cellinfo
    for i=1:cellnum1
        nth = i*2;
        F1cell(i).Position = [kC1(nth).Position(1)+v+5,kC1(nth).Position(2)+u+5];
        F1cell(i).String = kC1(nth).String;
        F1Contour(:,i) = reshape(circshift(kC1(nth-1).ZData,[u,v]),prod(SIZ),1);
    end
    %2nd image cellinfo
    for i=1:cellnum2
        nth = i*2;
        F2cell(i).Position = [kC2(nth).Position(1)-10,kC2(nth).Position(2)-10];
        F2cell(i).String = kC2(nth).String;
        F2Contour(:,i) = reshape(kC2(nth-1).ZData,prod(SIZ),1);
    end
else
     for ii=1:cellnum1
        nth = 2*cellnum1+1-ii*2;
        F1cell(ii).Position = [kC1(nth).Position(1)+v+5,kC1(nth).Position(2)+u+5];
        F1cell(ii).String = kC1(nth).String;
        tempmask= poly2mask(double(kC1(nth+1).XData),double(kC1(nth+1).YData),SIZ(1),SIZ(2));
        F1Contour(:,ii) = reshape(circshift(tempmask,[u,v]),[],1);
     end
     for ii=1:cellnum2
        nth = 2*cellnum2+1-ii*2;
        F2cell(ii).Position = [kC2(nth).Position(1)+5,kC2(nth).Position(2)+5];
        F2cell(ii).String = kC2(nth).String;
        tempmask2= poly2mask(double(kC2(nth+1).XData),double(kC2(nth+1).YData),SIZ(1),SIZ(2));
        F2Contour(:,ii) = reshape(tempmask2,[],1);
     end
end
%%
close all;

%%    %first image,red
fig3=figure()
imshowpair(fig1imgnew,fig2img)

title('fig1 in red & fig2 in green')

hold on;
for i=1:cellnum1
    %     text(F1cell(i).Position(1),F1cell(i).Position(2),F1cell(i).String,'Color',[ 1 0 0],'FontSize',8)
    contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,i),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[1 0 0])
end
%2nd image, green
for i=1:cellnum2
    %     text(F2cell(i).Position(1),F2cell(i).Position(2),F2cell(i).String,'Color',[ 0 1 0],'FontSize',8)
    contour(1:SIZ(2),1:SIZ(1),reshape(F2Contour(:,i),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[0 1 0])
end

%% Deselect unqualified cells first
x=[];y=[];
out1=zeros(1,cellnum1);
out2=zeros(1,cellnum2);
maskpic1 = zeros(SIZ);
maskpic2 = zeros(1,prod(SIZ));
% load(fullfile(path1a,[figname1(1:pos(end)) 'memmap.mat']),'magnification');
% Set potential good cells criteria: response and size
try 
    sbxread(fullfile(path1a,figname1(1:pos(3)-1)),0,1);
    global info;
    magnification = info.config.magnification;
    load('magnificationlist');
    mag= maglist(magnification)
catch
    disp('Update magnification otherwise mag=2'); 
    mag =2;
end
load('magnificationlist');

% cellsize = maglist(magnification)^2*30 *[.6 5] ; NOT USING CELLSIZE as a
% criteria at the moment JSUN Jan-2020
% cellsize = magnification^2*40 *4 ;
%%
while menu('Deselect?','Yes/Plot','Done')==1
    switch menu('How?','Rect','Pointers','Manuel input(w suggestion)','Plot')
        case 1
            rect = round(getrect);
            maskpic1(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3))=1;
            out1 = out1+reshape(maskpic1,1,prod(SIZ))*F1Contour;
            out2 = out2+reshape(maskpic1,1,prod(SIZ))*F2Contour;
        case 2
            [x,y]= ginput;
            maskpic2(sub2ind(SIZ, round(y),round(x)))=1;
            out1 = out1+maskpic2*F1Contour;
            out2 = out2+maskpic2*F2Contour;
        case 3
            prompt = {'1st red','2nd green'};
            total1 = F1Contour*(out1==0)';
            total2 = F2Contour*(out2==0)';
            %             total1=sum(F1Contour,2);
            %             total2=sum(F2Contour,2);
            defaultout1 = find(total2'*F1Contour==0); %| sum(F1Contour>.01)>cellsize(2)*4);
            defaultout2 = find(total1'*F2Contour==0); %| sum(F2Contour>.01)>cellsize(2)*4);
            default = {num2str(defaultout1),num2str(defaultout2)};
            figtitle = 'input bad cells';
            bad=inputdlg(prompt, figtitle,[2 80],default);
            out1(str2num(bad{1}))=1;
            out2(str2num(bad{2}))=1;
        case 4
            close(fig3)
            fig3 = figure()
             imshowpair(fig1imgnew,fig2img)
            set(gcf,'Position',[0 0 2000 1200])
            hold on;
            %first image,red
            for i=1:cellnum1
                if out1(i) ==0
                    text(F1cell(i).Position(1),F1cell(i).Position(2),F1cell(i).String,'Color',[ 1 0 0],'FontSize',8)
                    contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,i),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[1 0 0])
                end
            end
            %2nd image, green
            for i=1:cellnum2
                if out2(i) ==0
                    text(F2cell(i).Position(1),F2cell(i).Position(2),F2cell(i).String,'Color',[ 0 1 0],'FontSize',8)
                    contour(1:SIZ(2),1:SIZ(1),reshape(F2Contour(:,i),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[0 1 0])
                end
            end
    end
end
%
save([name '_cleaned.mat'],'out1','out2');
%% Select good cells and plot
prompt = 'red cell';
figtitle = 'input the green cell# matching the red';
pair1=[];
pair2=[];

figure()
imshowpair(fig1imgnew,fig2img)
set(gcf,'Position',[0 0 1000 600]);
hold on;

%2nd image, green
for i=1:cellnum2
    if out2(i) ==0
        text(F2cell(i).Position(1),F2cell(i).Position(2),F2cell(i).String,'Color',[ 0 1 0],'FontSize',8)
        contour(1:SIZ(2),1:SIZ(1),reshape(F2Contour(:,i),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[0 1 0])
    end
end
%%
i0=1; % the cell to start
% 1st image,red
low_overlap_temp = [];
while i0<=cellnum1
    if out1(i0) ==0
        text(F1cell(i0).Position(1),F1cell(i0).Position(2),F1cell(i0).String,'Color',[ 1 0 0],'FontSize',8)
        contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,i0),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[1 0 0])
        drawnow;
        default =find((F1Contour(:,i0)'*F2Contour) >0);
        [~,I] = sort(corr(F1Contour(:,i0),F2Contour(:,default)),'descend' ); %descending correlation number
        if corr(F1Contour(:,i0),F2Contour(:,default)) < 0.75
            low_overlap_temp(end+1) = i0;
        end
        general_overlap{i0} = corr(F1Contour(:,i0),F2Contour(:,default));
        %threshold for correlation to check after matching - display only
        %<75% for example
        gr=inputdlg({[prompt F1cell(i0).String],'Redo the last one'}, figtitle ,[2 80],{num2str(default(I)),'0'});
        if str2num(gr{2})==1
           i0 = i0-1;
        elseif   ~isempty(gr{1})
            i2=str2num(gr{1});
            temp= F1Contour(:,i0).*F2Contour(:,i2);
            if sum(temp(:)) == 0 || numel(i2)>1
                disp('wrong pairs? Redo matching'); i0 = i0-1;
            else
                pair1(end+1)=i0;
                pair2(end+1)=i2;
                text(F1cell(i0).Position(1),F1cell(i0).Position(2),F1cell(i0).String,'Color',[ 1 1 1],'FontSize',8)
                contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,i0),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[1 1 1])
                text(F2cell(i2).Position(1),F2cell(i2).Position(2),F2cell(i2).String,'Color',[ 1 1 1],'FontSize',8)
                contour(1:SIZ(2),1:SIZ(1),reshape(F2Contour(:,i2),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[1 1 1])
            end
        else
            text(F1cell(i0).Position(1),F1cell(i0).Position(2),F1cell(i0).String,'Color',[ 0 0 0],'FontSize',8)
            contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,i0),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[0 0 0])
        end
    end
        i0 = i0+1;
end

%% Modify selected pairs
while menu('ready to save?','Y','N')==2
     prompt = {'1st red','2nd green'};
     newpairs=inputdlg(prompt, 'make changes if desired',[2 80],{num2str(pair1),num2str(pair2)});
     pair1 = str2num(newpairs{1});
     pair2 = str2num(newpairs{2});
end
%% (redo step 3) Save all selected pairs
fig4=figure;
SIZ = size(fig1img);
fig1imgnew = circshift(fig1img,[u,v]);
imshowpair(fig1imgnew,fig2img)
set(gcf,'Position',[0 0 1200 1000]);
hold on;
%first image,red
for m=1:numel(pair1)
    text(F1cell(pair1(m)).Position(1),F1cell(pair1(m)).Position(2),F1cell(pair1(m)).String,'Color',[ 1 0 0],'FontSize',8)
    contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,pair1(m)),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[1 0 0])
    
    text(F2cell(pair2(m)).Position(1),F2cell(pair2(m)).Position(2),F2cell(pair2(m)).String,'Color',[ 0 1 0],'FontSize',8)
    contour(1:SIZ(2),1:SIZ(1),reshape(F2Contour(:,pair2(m)),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[0 1 0])
    
end

%%
selectedname=[name '_selected' num2str(numel(pair1))];
save(selectedname,'fig1img','fig2img','F1Contour','F2Contour','F1cell','F2cell','pair1','pair2','u','v');

saveas(fig4,selectedname,'fig');
saveas(fig4,selectedname,'png');

close(fig4)
% extra_name =[name '_additional_information.mat'];
% save(extra_name, 'green_duplicates', 'red_indices', 'general_overlap');
%% plot &compare twodays
% plot_twoday_data
