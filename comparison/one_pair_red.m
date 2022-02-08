function one_pair_red(a,b,c)
%a = red cells 
%b = cells
%c = all pairs pair
%one_pair_red(find(se1{1}.redcell(:,1)==1),find(se1{1}.iscell(:,1)==1),allpairs{1}.pair1')
%one_pair_red(find(se1{1}.redcell2(:,4)==1),find(se1{1}.iscell(:,1)==1),allpairs{1}.pair1')
    
% b = find(se1{ii}.iscell(:,1)==1);
%    
% redinpairidx_1 = find(se1{ii}.redcell2(iscellidx_1(allpairs{ii}.pair1))==1);
% sumredinpair1 = sum(se1{ii}.redcell2(iscellidx_1(allpairs{ii}.pair1)))
% 
num_red_cells = 0;
for i = 1:length(a)
    for p = 1:length(b)
        if a(i) == b(p)
            for j = 1:length(c)
                if c(j) == p
                    num_red_cells = num_red_cells + 1;
%                     fprintf('i, m(i), l(p), o(j)')
%                     i
%                     m(i)
%                     l(p)
%                     p
%                     j
                else
                    continue
                end
            end
        end
    end
end
num_red_cells
end