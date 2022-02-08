function img=deletetxt(img)
        axesObjs = get(img, 'Children');  %axes handles
        dataObjs = get(axesObjs, 'Children');
for i=1:size(dataObjs)
    if strcmp(dataObjs(i).Type,'text')
        delete(dataObjs(i));
    end
end
%figure;plot(dataObjs(101).ContourMatrix(1,2:end),dataObjs(101).ContourMatrix(2,2:end))