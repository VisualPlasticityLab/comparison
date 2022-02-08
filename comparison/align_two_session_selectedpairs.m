%%
fn=dir('*selected*.mat');
load(fn.name);
%%
figure;
SIZ = size(fig1img);
fig1imgnew = circshift(fig1img,[u,v]);
imshowpair(fig1imgnew,fig2img)
set(gcf,'Position',[0 0 1200 1000]);
hold on;
%first image,red
selectedpairs= 1:10;

for m=selectedpairs
    text(F1cell(pair1(m)).Position(1),F1cell(pair1(m)).Position(2),F1cell(pair1(m)).String,'Color',[ 1 0 0],'FontSize',8)
    contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,pair1(m)),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[1 0 0])
    
    text(F2cell(pair2(m)).Position(1),F2cell(pair2(m)).Position(2),F2cell(pair2(m)).String,'Color',[ 0 1 0],'FontSize',8)
    contour(1:SIZ(2),1:SIZ(1),reshape(F2Contour(:,pair2(m)),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[0 1 0])
    
end