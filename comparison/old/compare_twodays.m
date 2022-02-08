%function compare_twodays

%% load the files and figs
[figname1,path1a]=uigetfile('.fig','find the 1st fig');
% [~,path1b]=uigetfile('peakSI.mat','locate the 1st peakSI.mat file');
% A=load(fullfile(path1a,'peakSI.mat'));
[figname2,path2a]=uigetfile('.fig','find the 2nd fig');
% [~,path2b]=uigetfile('peakSI.mat','locate the 2nd peakSI.mat file');
% B=load(fullfile(path2a,'peakSI.mat'));
openfig(fullfile(path1a,figname1));
fig1 = gca;
fig1img =getimage(fig1);
kC= fig1.Children;
cellnum=(numel(kC)-1)/2


openfig(fullfile(path2a,figname2));
fig2 = gca;
fig2img= getimage(fig2);
% fig1img=img;
%fig1img=img(:,:,2);
%%
fig1img = fig1img(:,:,1);
fig2img = fig2img(:,:,1);

fig1img=fig1img*(max(fig2img(:))/max(fig1img(:)));
[u v] = fftalign(fig1img,fig2img);
fig1imgnew = circshift(fig1img,[u,v]);
% figure;imshowpair(fig1img,fig2img);
% [optimizer, metric] = imregconfig('multimodal');
% [fig1img_reg,R_reg] =  imregister(fig1img, fig2img, 'translation', optimizer, metric);
% tform = imregtform(fig1img, fig2img, 'translation', optimizer, metric);
%fig1imgRegistered= imwarp(fig1img,tform,'OutputView',imref2d(size(fig1img)));

figure;
imshowpair(fig1imgnew,fig2img)
for i=1:cellnum
    text(kC(i).Position(1)+u,kC(i).Position(2)+v,kC(i).String,'Color',[ 0 0 1],'font',10)
end
for j=cellnum+1:cellnum*2
    contour(kC(j).XData,kC(j).YData,circshift(kC(j).ZData,[u v]),[0.01 1])
end

name=[strtok(figname1(1:end-4),'cell') '&' strtok(figname2(1:end-4),'cell')];
title(strrep(name,'_','-'));
% saveas(fig2,name,'fig');
% saveas(fig2,name,'png');
%% align the figures and find correponding cells
% saveas(fig1,fullfile(path1a,[strtok(figname1,'.') '.png']),'png');
% fig1=deletetxt(fig1);
% figname1_notxt=fullfile(path1a,[strtok(figname1,'.') 'notxt.png']);
% saveas(fig1,figname1_notxt,'png');
% 
% fig1=deletelabel(fig1);
% figname1_nolabel=fullfile(path1a,[strtok(figname1,'.') 'nolabel.png']);
% saveas(fig1,figname1_nolabel,'png');
% 
% saveas(fig2,fullfile(path2a,[strtok(figname2,'.') '.png']),'png');
% fig2=deletetxt(fig2);
% figname2_notxt=fullfile(path2a,[strtok(figname2,'.') 'notxt.png']);
% saveas(fig2,figname2_notxt,'png');
% 
% fig2=deletelabel(fig2);
% figname2_nolabel=fullfile(path2a,[strtok(figname2,'.') 'nolabel.png']);
% saveas(fig2,figname2_nolabel,'png');

moving=getimage(fig1);
fixed=getimage(fig2);
moving=moving*(max(fixed(:))/max(moving(:)))/10000;
fixed = fixed/10000;

figure;imshowpair(moving,fixed);
[optimizer, metric] = imregconfig('multimodal');
[moving_reg,R_reg] =  imregister(moving, fixed, 'affine', optimizer, metric);
tform = imregtform(moving, fixed, 'affine', optimizer, metric);
disp(tform.T)
movingRegistered= imwarp(moving,tform,'OutputView',imref2d(size(moving)));%have to set 'OutputView' mode
figure;imshowpair(movingRegistered,fixed);

%defaultpair=num2cell([long,out]);
selected=inputdlg('input matched cell#','Exclusion',[4 80],defaultpair);
pair=str2num(selected{1});

%%


%% plot the figures and save data for future analysis?

for i=1:size(pair,1)
    sig(:,:,:,2*i-1)=A.SI.sigF(:,:,:,find(A.Cor==pair(i,1)));
    
    sig(:,:,:,2*i)=B.SI.sigF(:,:,:,find(B.Cor==pair(i,2)));
end
Cor=reshape(pair',2,[]);
win_sig=A.win_sig;
matrix=[];
hsigF=sigFplt(sig,matrix,win_sig,Cor);  % sigF:seg,rep,Var,ncell
