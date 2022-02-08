function countred(se1,se2,allpairs)

fprintf('redchannel correction by defult\n')

for ii=1:3 
    iscellidx_1 = find(se1{ii}.iscell(:,1)==1)'; 
    %find only the cells out of all ROI
    prered=sum(se1{ii}.redcell(iscellidx_1(allpairs{ii}.pair1),1)==1);
    iscellidx_2 = find(se2{ii}.iscell(:,1)==1)'; 
    postred = sum(se2{ii}.redcell(iscellidx_2(allpairs{ii}.pair2),1)==1);
    bothred = sum(se1{ii}.redcell(iscellidx_1(allpairs{ii}.pair1),1)==1&se2{ii}.redcell(iscellidx_2(allpairs{ii}.pair2),1)==1);

    fprintf('plane%d:of%dpairs, red: pre-stim %d,post-stim %d,both %d\n',...
        ii-1,numel(allpairs{ii}.pair1),prered,postred,bothred)

end


fprintf('redchannel correction by JS\n')

for ii=1:3 
    iscellidx_1 = find(se1{ii}.iscell(:,1)==1)'; 
    iscellidx_2 = find(se2{ii}.iscell(:,1)==1)';
    %find only the cells out of all ROI
    prered=sum(se1{ii}.redcell2(iscellidx_1(allpairs{ii}.pair1),4)==1);
    postred = sum(se2{ii}.redcell2(iscellidx_2(allpairs{ii}.pair2),4)==1);
    bothred = sum(se1{ii}.redcell2(iscellidx_1(allpairs{ii}.pair1),4)==1&se2{ii}.redcell2(iscellidx_2(allpairs{ii}.pair2),4)==1);

    fprintf('plane%d:of%dpairs, red: pre-stim %d,post-stim %d,both %d\n',...
        ii-1,numel(allpairs{ii}.pair1),prered,postred,bothred)

end