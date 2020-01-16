function boundaryNodeIDs = getGlobalBoundaryNodeIDs(eCONN,eDEGREE)
boundaryNodeIDs = [];
for e = 1:size(eCONN,2)
    elemBoundaryMap = buildElementBoundaryMap(eDEGREE(:,e));
    for s = 1:length(elemBoundaryMap)
        nSharedNodes = sum(ismember(eCONN,eCONN(elemBoundaryMap{s},e)));
        nSharedFaces = sum(nSharedNodes == length(elemBoundaryMap{s}))-1;
        if nSharedFaces == 0
            boundaryNodeIDs = [boundaryNodeIDs; eCONN(elemBoundaryMap{s},e)];
        end
    end
end