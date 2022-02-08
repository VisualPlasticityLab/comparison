%%compare two analyzed image from memmapfiles

for i=1:2
    [f,p]=uigetfile('.signals','load calcium signal data .signals');
    f=fullfile(p,f);
    data=loadjson( [f(1:findstr(f,'cell')+3) '.jmesh']);
    F{i}.circle= data.jmesh{:};
    fn=f(1:findstr(f,'_memmap')-1);
    F{i}.m=load([fn '.align'],'m','-mat');
    F{i}.T=load([fn '.align'],'T','-mat');
end

for i=1:numel(circle)
    cntr(i,:)=circle{i}.centroid;
end

for i=1:numel(circle_1)
    cntr1(i,:)=circle_1{i}.centroid;
end

load('240_915_001.align','-mat')
m1=m;T1=T;
load('240_915_002.align','-mat')
m2=m;T2=T;
[u v] = fftalign(m1,m2);
imshowpair(m1,m2)


figure;hold on;
scatter(cntr(:,2),cntr(:,1),'*');
text(cntr(:,2),cntr(:,1),num2str([1:numel(circle)]'),'Color','b')
scatter(cntr1(:,2),cntr1(:,1),'r*')
text(cntr1(:,2),cntr1(:,1),num2str([1:numel(circle_1)]'),'Color','r')
axis ij