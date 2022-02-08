% function align_images
% this function compare two images and find physically matched cell, check
% cell responses with the function analysis_twoday_data, see more in
% Demo_2p
close all;
clear all;clc
%% load two files and figs
[aligned1,path1a]=uigetfile('.align','find the 1st alignment');
[figname1,path1a]=uigetfile([path1a '*.fig'],'find the 1st figname');
m1 = load(fullfile(path1a,aligned1),'-mat');
pos1=strfind(strtok(figname1,'c'),'_');
figplane1 = str2num(figname1(pos1(end)-1));

[aligned2,path2a]=uigetfile('.align','find the 2nd alignment');
[figname2,path2a]=uigetfile([path2a '.fig'],'find the 2nd figname');
m2 = load(fullfile(path2a,aligned2),'-mat');
pos2= strfind(strtok(figname2,'c'),'_');
figplane2 = str2num(figname2(pos2(end)-1));
name=[figname1(1:pos1(end)-1) '&' figname2(1:pos2(end)-1) ];
%% align two images and plot figure 1 in green & figure2 in purple
if isfield(m1,'mg') & isfield(m2,'mg')
    align1img= squeeze(m1.mg(:,:,figplane1));
    align2img= squeeze(m2.mg(:,:,figplane2));
else
    align1img= squeeze(m1.m(:,:,figplane1));
    align2img= squeeze(m2.m(:,:,figplane2));
end

fig3 = figure();
while menu('Select areas for image alignment?','Yes','Done')==1
    imshowpair(align1img,align2img);
    rect=round(getrect(fig3));
    img1=align1img(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3));
    img2=align2img(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3));
    [u v] = fftalign(img1,img2);
    align1imgnew = circshift(align1img,[u,v]);
    imshowpair(align1imgnew,align2img);
    imgcorr = corr(double(align1imgnew(:)),double(align2img(:)));
    img1new = align1imgnew(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3));
    imgcorr2 = corr(double(img1new(:)),double(img2(:)));
    title(sprintf('u=%d,v=%d,corr=%.2f,select corr=%.2f',u,v,imgcorr,imgcorr2))
    %     set(fig3,'Position',[0 0 2000 1200]);
end
%% 1st image cellinfo
try
    fig1=openfig(fullfile(path1a,figname1));
    fig1img =getimage(fig1);
catch
    [day1sig,p1]=uigetfile('.signals','select signal file for day1eye');
    s1 = matfile(fullfile(p1,day1sig));
    fig1img =s1.V;
    %s1.A has all the contours and text position is   com(Aor(:,1:end),d1,d2);
    fig1=figure();
    plot_contours(s1.A,fig1img,0.95,1); %thr=.95,display=TRUE
    saveas(fig1,[figname1 '.fig'],'fig');
    saveas(fig1,[figname1 '.png'],'png');
end
kC1= fig1.Children.Children;
cellnum1=(numel(kC1)-1)/2;

SIZ = size(fig1img);
F1Contour =zeros(prod(SIZ),cellnum1);
for i=1:cellnum1
    nth = cellnum1+1-i;
    F1cell(i).Position = [kC1(nth).Position(1)+v+5,kC1(nth).Position(2)+u+5];
    F1cell(i).String = kC1(nth).String;
    F1Contour(:,i) = reshape(circshift(kC1(nth+cellnum1).ZData,[u,v]),prod(SIZ),1);
end
fig1imgnew = circshift(fig1img,[u,v]);
%% 2nd image cellinfo
try
    fig2=openfig(fullfile(path2a,figname2));
    fig2img= getimage(fig2);
catch
    [day2sig,p2]=uigetfile('.signals','select signal file for day2eye');
    s2 = matfile(fullfile(p2,day2sig));
    fig2img =s2.V;
    fig2=figure();
    plot_contours(s2.A,s2.V,0.95,1);
    saveas(fig2,[figname2 '.fig'],'fig');
    saveas(fig2,[figname2 '.png'],'png');
end
kC2= fig2.Children.Children;
cellnum2=(numel(kC2)-1)/2;
F2Contour =zeros(prod(SIZ),cellnum2);
for i=1:cellnum2
    nth = cellnum2+1-i;
    F2cell(i).Position = [kC2(nth).Position(1)-10,kC2(nth).Position(2)-10];
    F2cell(i).String = kC2(nth).String;
    F2Contour(:,i) = reshape(kC2(nth+cellnum2).ZData,prod(SIZ),1);
end
save([name  '_aligned.mat'],'name','cellnum1','cellnum2','u','v','fig1imgnew','fig2img','F1cell','F2cell','F1Contour','F2Contour','-v7.3');
%% Select GD cells from each days first
try
    load(fullfile(path1a,[figname1(1:pos1(end)) 'memmap.mat']),'magnification');
    if mod(magnification,1)>0 || magnification ==2
        magnification = magnification*2.5;
    end
    load magnificationlist.mat ;
    mag = maglist(magnification);
catch
    mag = 2;
end
[f1.pair0,f1.pair1,f1.pair2,c1]=autoselect(fig1imgnew,F1Contour,mag);%f1.pair0 selected area, f1.pair1 selected shape&response
[f2.pair0,f2.pair1,f2.pair2,c2]=autoselect(fig2img,F2Contour,mag);
f1.selected = find(f1.pair0&f1.pair1&f1.pair2);
f2.selected = find(f2.pair0&f2.pair1&f2.pair2);

fig3 = figure();
imshowpair(fig1imgnew,fig2img) %first image in green and second image in purple
hold on
for i=f1.selected
    contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,i),SIZ(1),SIZ(2)),[0.01 1],...
        'LineColor',[1 0 0])                    
    text(F1cell(i).Position(1),F1cell(i).Position(2),F1cell(i).String,'Color',[ 1 0 0],'FontSize',8)
end
for i=f2.selected
    contour(1:SIZ(2),1:SIZ(1),reshape(F2Contour(:,i),SIZ(1),SIZ(2)),[0.01 1],...
        'LineColor',[0 0 1],'Linewidth',0.5)
    text(F2cell(i).Position(1),F2cell(i).Position(2),F2cell(i).String,'Color',[ 0 0 1],'FontSize',8)
end
saveas(fig3,sprintf('%s_%d&%d_aligned.fig',...
    name,numel(f1.selected),numel(f2.selected)))
%%
x=[];y=[];
cellnum1 = size(F1Contour,2);
cellnum2 = size(F2Contour,2);
out1=~(f1.pair0&f1.pair1);
out2=~(f2.pair0&f2.pair1);
maskpic1 = zeros(SIZ);
maskpic2 = zeros(1,prod(SIZ));

while menu('Deselect?','Yes/Plot','Done')==1
    switch menu('How?','Rect','Pointers','Manuel input(w suggestion)','Plot','Save')
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
            total1 = F1Contour*(out1==0)';%overall area with img1 cells
            total2 = F2Contour*(out2==0)';%overall area with img2 cells
           
            % select rule (stricter): cells not overlapping more than 50% with any cells in the
            % other image
            defaultout1 = find((double(total2'>0)*double(F1Contour>.01)./sum(F1Contour>.01,1))<.5 | sum(F1Contour>.01,1)>cellsize*3);
            defaultout2 = find((double(total1'>0)*double(F2Contour>.01)./sum(F2Contour>.01,1))<.5 | sum(F2Contour>.01)>cellsize*3);
            
            default = {num2str(defaultout1),num2str(defaultout2)};
            prmtitle = 'input bad cells';
            bad=inputdlg(prompt, prmtitle,[2 80],default);
            out1(str2num(bad{1}))=1;
            out2(str2num(bad{2}))=1;
        case 4
            clf;
            imshowpair(fig1imgnew,fig2img)
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
                    text(F2cell(i).Position(1),F2cell(i).Position(2),F2cell(i).String,'Color',[ 0 0 1],'FontSize',8)
                    contour(1:SIZ(2),1:SIZ(1),reshape(F2Contour(:,i),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[0 0 1])
                end
            end
        case 5
%             saveas(gcf,[name '_cleaned'],'fig');
%             saveas(gcf,[name '_cleaned.png'],'png');
            save([name '_cleaned.mat'],'out1','out2');
    end
end
close all
%% Select good cells and plot
pair1=[];
pair2=[];
prompt = 'red cell';
prmtitle = 'input the green cell# matching the red';
figure
imshowpair(fig1imgnew,fig2img)
set(gcf,'Position',[0 0 2000 1200]);
hold on;

i=1;

while i<=cellnum1
    if out1(i) ==0
        text(F1cell(i).Position(1),F1cell(i).Position(2),F1cell(i).String,'Color',[ 1 0 0],'FontSize',8)
        contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,i),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[1 0 0])
        drawnow;
        default =find((F1Contour(:,i)'*F2Contour) >0 &~out2);
        [~,I] = sort(corr(F1Contour(:,i),F2Contour(:,default)),'descend' );
        if numel(default(I))
            for k=1:numel(default(I))
                text(F2cell(default(k)).Position(1),F2cell(default(k)).Position(2),F2cell(default(k)).String,'Color',[ 0 0 1],'FontSize',8)
                contour(1:SIZ(2),1:SIZ(1),reshape(F2Contour(:,default(k)),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[0 0 1])
            end
        end
        gr=inputdlg({[prompt num2str(i)]}, prmtitle ,[2 80],{num2str(default(I))});
        if ~isempty(gr{1})
            i2=str2num(gr{1});
            temp= F1Contour(:,i).*F2Contour(:,i2);
            if sum(temp(:)) == 0
                disp('wrong pairs? Redo matching'); i = i-1;
            elseif numel(i2)>1
                disp('cannot have one than more match');i = i-1;
            else
                pair1(end+1)=i;
                pair2(end+1)=i2;
            end
        end
        text(F1cell(i).Position(1),F1cell(i).Position(2),F1cell(i).String,'Color',[ 1 1 1],'FontSize',8)
        contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,i),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[1 1 1])
        
        
        for k=1:numel(default(I))
            text(F2cell(default(k)).Position(1),F2cell(default(k)).Position(2),F2cell(default(k)).String,'Color',[ 1 1 1],'FontSize',8)
            contour(1:SIZ(2),1:SIZ(1),reshape(F2Contour(:,default(k)),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[1 1 1])
        end
    end
    i = i+1;
end
close(gcf)

%%
fig4=figure('Position',[0 0 2000 1200]);
imshow(overlayfigs(fig1imgnew,fig2img))
hold on;
%first cells,red
for m=1:numel(pair1)
    text(F1cell(pair1(m)).Position(1),F1cell(pair1(m)).Position(2),F1cell(pair1(m)).String,'Color',[ 1 0 0],'FontSize',8)
    contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,pair1(m)),SIZ(1),SIZ(2)),[0.01 1],...
        'LineColor',[1 0 0],'LineWidth',.5)
    
     text(F2cell(pair2(m)).Position(1),F2cell(pair2(m)).Position(2),F2cell(pair2(m)).String,'Color',[ 0 0 1],'FontSize',8)
    contour(1:SIZ(2),1:SIZ(1),reshape(F2Contour(:,pair2(m)),SIZ(1),SIZ(2)),[0.01 1],...
        'LineColor',[0 0 1],'LineWidth',.5)
end
%%
selectedname=sprintf('%s_%d&%d',name,numel(f1.selected),numel(f2.selected))
copyfile([name  '_aligned.mat'],[selectedname '.mat']);

pair1=[];
pair2=[];
for i=f1.selected
    default =find(F1Contour(:,i)'*F2Contour(:,f2.selected) >0 );
    if ~isempty(default)
        prompt = sprintf('image1 cell#%d : image2 cell#%d,confirm?',i,f2.selected(default));
        if strcmp(input(prompt,'s'),'y')
        pair1(end+1) = i;pair2(end+1) = f2.selected(default);
        end
    end
end
save([selectedname '.mat'],'SIZ','mag','pair1','pair2','f1','f2','-append');
%%
hh = figure('Position',[200 400 452 380]);
subplot(1,3,1)
imagesc(fig1imgnew);axis off
colormap gray; hold on
for i=find(f1.pair0&f1.pair1)
    contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,i),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[1 0 0])
end

subplot(1,3,2)
imagesc(fig2img);axis off
colormap gray; hold on
for i=find(f2.pair0&f2.pair1)
    contour(1:SIZ(2),1:SIZ(1),reshape(F2Contour(:,i),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[0 0 1])
end

subplot(1,3,3)
imagesc(overlayfigs(fig1imgnew,fig2img));axis off
hold on;
for i=pair1
    contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,i),SIZ(1),SIZ(2),[]),[0.01 1],'LineColor',[1 0 0])
    text(F1cell(i).Position(1),F1cell(i).Position(2),F1cell(i).String,'Color',[ 1 0 0],'FontSize',8)
end
%2nd image, Blue
for i=pair2
    contour(1:SIZ(2),1:SIZ(1),reshape(F2Contour(:,i),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[0 0 1])
    text(F2cell(i).Position(1),F2cell(i).Position(2),F2cell(i).String,'Color',[0 0 1],'FontSize',8)
end
%%
saveas(hh,[selectedname '.png']);
saveas(hh,[selectedname '.fig']);

%% optional 
% plotODIcomparison
% analysis_twoday_data
% plotposition
% cacorrelation
