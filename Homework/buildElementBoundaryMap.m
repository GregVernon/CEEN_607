function SS = buildElementBoundaryMap(eDegree)
nElemNodes = eDegree + 1;
nDimensions = length(eDegree);

if nDimensions == 1
    SS = cell(2,1);
    SS{1} = 1;
    SS{2} = nElemNodes;
elseif nDimensions == 2
    SS = cell(4,1);
    NX = nElemNodes(1);
    NY = nElemNodes(2);
    nodeID = 0;
    for nj = 1:NY
        for ni = 1:NX
            nodeID = nodeID + 1;
            if ni == 1
                SS{1} = [SS{1}; nodeID];
            end
            if ni == NX
                SS{2} = [SS{2}; nodeID];
            end
            if nj == 1
                SS{3} = [SS{3}; nodeID];
            end
            if nj == NY
                SS{4} = [SS{4}; nodeID];
            end
        end
    end
    
elseif nDimensions == 3
    SS = cell(6,1);
    NX = nElemNodes(1);
    NY = nElemNodes(2);
    NZ = nElemNodes(3);
    nodeID = 0;
    for nk = 1:NZ
        for nj = 1:NY
            for ni = 1:NX
                nodeID = nodeID + 1;
                if ni == 1
                    SS{1} = [SS{1}; nodeID];
                end
                if ni == NX
                    SS{2} = [SS{2}; nodeID];
                end
                if nj == 1
                    SS{3} = [SS{3}; nodeID];
                end
                if nj == NY
                    SS{4} = [SS{4}; nodeID];
                end
                if nk == 1
                    SS{5} = [SS{5}; nodeID];
                end
                if nk == NZ
                    SS{6} = [SS{6}; nodeID];
                end
            end
        end
    end
end
end