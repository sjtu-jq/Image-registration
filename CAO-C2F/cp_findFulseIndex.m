function fulseindex = cp_findFulseIndex(correctindex,length)
% Find fulse indices of the candidate matches
fulseindex = [];
for i=1:length
    if isempty(find(correctindex == i))
        fulseindex = [fulseindex;i];
    end
end

