% function align_images
% this function compare two images and find physically matched cell, check
% cell responses with the function analysis_twoday_data, see more in
% Demo_2p
close all;
clear all
%% load the files and figs

[aligned1,path1a]=uigetfile('.align','find the 1st alignment');
[aligned2,path2a]=uigetfile('.align','find the 2nd alignment');

m1 = load(fullfile(path1a,aligned1),'-mat');
m2 = load(fullfile(path2a,aligned2),'-mat');

pos_ = strfind(aligned1,'_');
day1 = aligned1(pos_(1)+1:pos_(2)-1);
day2 = aligned2(pos_(1)+1:pos_(2)-1);
name=[day1 '&' day2];
%% align two images and plot figure 1 in green & figure2 in purple
if isfield(m1,'mg') & isfield(m2,'mg')
    align1img= mean(m1.mg,3);
    align2img= mean(m2.mg,3);
else
    align1img= mean(m1.m,3);
    align2img= mean(m2.m,3);
end

[u(1,1) v(1,1)]=fftalignimg(m1.m(:,:,1),m1.m(:,:,2));
[u(1,3) v(1,3)]=fftalignimg(m1.m(:,:,3),m1.m(:,:,2));

[u(2,1) v(2,1)]=fftalignimg(m2.m(:,:,1),m1.m(:,:,1));
[u(2,2) v(2,2)]=fftalignimg(m2.m(:,:,2),m1.m(:,:,2));
[u(2,3) v(2,3)]=fftalignimg(m2.m(:,:,3),m1.m(:,:,3));

for pln=1:3
img1(:,:,pln)=circshift(m1.m(:,:,pln),[u(1,pln) v(1,pln)]);
img2(:,:,pln)=circshift(m2.m(:,:,pln),[u(1,pln) v(1,pln)]+[u(2,pln) v(2,pln)]);
end

img1 = double(img1);
img2 = double(img2);

fig1img = img1/max(reshape(img1(200:end-200,200:end-200,:),[],1));
fig2img = img2/max(reshape(img2(200:end-200,200:end-200,:),[],1));
fig3 = figure();imshow(fig1img.*(1-fig2img))

%% Take cell information
day1file = matfile([day1 '\allplanes.mat']);
day2file = matfile([day2 '\allplanes.mat']);

cellnum1 = day1file.cellnum;
level_1 = discretize(1:sum(cellnum1),[1 cumsum(cellnum1)+1]);
SIZ = size(fig1img);
tempcell = day1file.F1cell;
tmpcontour = reshape(day1file.F1Contour,SIZ(1),SIZ(2),[]);
F1cell.String = tempcell.String;

for ii=1:sum(cellnum1)
    F1cell(ii).Position = tempcell(ii).Position +[u(1,level_1(ii)) v(1,level_1(ii))];
    F1Contour(:,ii) = reshape(circshift(day1file.F1Contour(:,ii),[u(1,level_1(ii)) v(1,level_1(ii))],prod(SIZ),1);
end
%%
saveas(fig3,[name  '_aligned.png'])
save([name  '_aligned'],'u','v','fig1img','fig2img','F1cell','F2cell','F1Contour','F2Contour','-v7.3');
%% first image,green
fig3 = figure('Position',[0 0 2000 1200]);
SIZ = size(fig2img);
imshowpair(fig1imgnew,fig2img)
hold on;
for i=1:numel(F1cell)
    %     text(F1cell(i).Position(1),F1cell(i).Position(2),F1cell(i).String,'Color',[ 1 0 0],'FontSize',8)
    contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,i),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[1 0 0])
end

%2nd image, green
for i=1:numel(F2cell)
    %     text(F2cell(i).Position(1),F2cell(i).Position(2),F2cell(i).String,'Color',[ 0 1 0],'FontSize',8)
    contour(1:SIZ(2),1:SIZ(1),reshape(F2Contour(:,i),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[0 1 0])
end

% saveas(fig3,[name '_overlapped'],'fig');
% saveas(fig3,[name '_overlapped'],'png');
%% Deselect unqualified cells first
try
    load(fullfile(path1a,[figname1(1:pos(end)) 'memmap.mat']),'magnification');
    if mod(magnification,1)>0 || magnification ==2
        magnification = magnification*2.5;
    end
catch
    magnification =5;
end
%cellsize = maglist(magnification*2.5)* 40 *4 ;

load magnificationlist.mat ;
cellsize = maglist(magnification)^2* 60  ;%%

%%
x=[];y=[];
out1=zeros(1,cellnum1);
out2=zeros(1,cellnum2);
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
                    text(F2cell(i).Position(1),F2cell(i).Position(2),F2cell(i).String,'Color',[ 0 1 0],'FontSize',8)
                    contour(1:SIZ(2),1:SIZ(1),reshape(F2Contour(:,i),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[0 1 0])
                end
            end
        case 5
%             saveas(gcf,[name '_cleaned'],'fig');
            saveas(gcf,[name '_cleaned.png'],'png');
            save([name '_cleaned.mat'],'out1','out2');
    end
end
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
%%
while i<=cellnum1
    if out1(i) ==0
        text(F1cell(i).Position(1),F1cell(i).Position(2),F1cell(i).String,'Color',[ 1 0 0],'FontSize',8)
        contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,i),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[1 0 0])
        drawnow;
        default =find((F1Contour(:,i)'*F2Contour) >0 &~out2);
        [~,I] = sort(corr(F1Contour(:,i),F2Contour(:,default)),'descend' );
        if numel(default(I))
            for k=1:numel(default(I))
                text(F2cell(default(k)).Position(1),F2cell(default(k)).Position(2),F2cell(default(k)).String,'Color',[ 0 1 0],'FontSize',8)
                contour(1:SIZ(2),1:SIZ(1),reshape(F2Contour(:,default(k)),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[0 1 0])
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

fig4=figure(4);
imshowpair(fig1imgnew,fig2img)
set(gcf,'Position',[0 0 2000 1200]);
hold on;
%first image,red
for m=1:numel(pair1)
    text(F1cell(pair1(m)).Position(1),F1cell(pair1(m)).Position(2),F1cell(pair1(m)).String,'Color',[ 1 0 0],'FontSize',8)
    contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,pair1(m)),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[1 0 0])
    
    text(F2cell(pair2(m)).Position(1),F2cell(pair2(m)).Position(2),F2cell(pair2(m)).String,'Color',[ 0 1 0],'FontSize',8)
    contour(1:SIZ(2),1:SIZ(1),reshape(F2Contour(:,pair2(m)),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[0 1 0])
    
end
%%
hh = figure('Position',[0 0 2000 1200]);
imshowpair(fig1imgnew,fig2img)
set(gcf,'Position',[0 0 2000 1200]);
hold on;
for i=pair1
    contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,i),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[1 0 0])
end
%2nd image, green
for i=pair2
    contour(1:SIZ(2),1:SIZ(1),reshape(F2Contour(:,i),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[0 1 0])
end

selectedname=[name '_selected' num2str(numel(pair1))]
save(selectedname,'pair1','pair2');
saveas(hh,[selectedname '.png']);

%% optional 
% plotODIcomparison
% analysis_twoday_data
% plotposition
% cacorrelation
