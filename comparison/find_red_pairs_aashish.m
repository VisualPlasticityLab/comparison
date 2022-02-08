%aashish
%function [sumredinpair1, sumredinpair2, num_red_matched] = find_red_pairs(se1, se2, allpairs)
function find_red_pairs_aashish(se1, se2, allpairs)
%finds red matched pairs 
% gets first plane right, others no
for ii = 1:3
    %ii = plane number
    ii;
    %prestim
    iscellidx_1 = find(se1{ii}.iscell(:,1)==1);
    iscell_1(ii) = sum(se1{ii}.iscell(:,1)==1);
    
    %poststim
    iscellidx_2 = find(se2{ii}.iscell(:,1)==1);
    iscell_2(ii) = sum(se2{ii}.iscell(:,1)==1);
    
%     %uses redcell (suite2p)
    channel = ' '
    redinpairidx_1 = find(se1{ii}.redcell(iscellidx_1(allpairs{ii}.pair1))==1);
    sumredinpair1(ii) = sum(se1{ii}.redcell(iscellidx_1(allpairs{ii}.pair1)));

    redinpairidx_2 = find(se2{ii}.redcell(iscellidx_2(allpairs{ii}.pair2))==1);
    sumredinpair2(ii) = sum(se2{ii}.redcell(iscellidx_2(allpairs{ii}.pair2)));
    
% %     %uses redcell2 (jennifer)
%     channel = '2';
%     redinpairidx_1 = find(se1{ii}.redcell2(iscellidx_1(allpairs{ii}.pair1),4)==1);
%     sumredinpair1(ii) = sum(se1{ii}.redcell2(iscellidx_1(allpairs{ii}.pair1),4));
% 
%     redinpairidx_2 = find(se2{ii}.redcell2(iscellidx_2(allpairs{ii}.pair2),4)==1);
%     sumredinpair2(ii) = sum(se2{ii}.redcell2(iscellidx_2(allpairs{ii}.pair2),4));
   
    matches = [];
    for num_1 = 1:length(redinpairidx_1)
        %num_1
        for num_2 = 1:length(redinpairidx_2)
           % num_2
            if redinpairidx_1(num_1) == redinpairidx_2(num_2)
                %redinpairidx_1
                %redinpairidx_2
                matches(end+1) = redinpairidx_1(num_1);
            end
        end
    end
    
    matches  %array of red matched cell index numbers
    num_red_matched(ii) = length(matches); %number of matched pairs
   	%fprintf('correction: redcell%cchannel') fix
    fprintf('\nplane%d:of %dpairs, cells: prestim %d, poststim %d, red: pre-stim %d,post-stim %d,both %d\n',...
        ii-1,numel(allpairs{ii}.pair1),iscell_1(ii), iscell_2(ii), sumredinpair1(ii),sumredinpair2(ii),num_red_matched(ii))
    
end

end