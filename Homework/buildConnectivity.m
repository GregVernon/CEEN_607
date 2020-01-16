function [eCONN,eDegree] = buildConnectivity(nElemNodes)
if iscell(nElemNodes) == true
    nDimensions = length(nElemNodes);
    nElem = cellfun(@length,nElemNodes);
else
    nDimensions = 1;
    nElem = length(nElemNodes);
end

eDegree = zeros(nDimensions,prod(nElem));
if nDimensions == 1
    eCONN = nan(max(nElemNodes),nElem);
    for e = 1:nElem
        eDegree(1,e) = nElemNodes(e)-1;
        if e == 1
            nodeIDs = 1:nElemNodes(e);
        else
            nodeIDs = nodeIDs(end):(nodeIDs(end)+(nElemNodes(e)-1));
        end
        eCONN(1:nElemNodes(e),e) = nodeIDs;
    end
elseif nDimensions == 2
    nTotElem = prod(nElem);
    maxElemNodes = max(max(kron(nElemNodes{1},nElemNodes{2}')));
    eCONN = nan(maxElemNodes,nTotElem);
    
    nxNodes = (nElem(1)+1) + sum(nElemNodes{1}-2);
    nyNodes = (nElem(2)+1) + sum(nElemNodes{2}-2);
    
    NX = nElem(1);
    NY = nElem(2);
    for ey = 1:NY
        for ex = 1:NX
            e = (ey-1)*NY + ex;
            eDegree(:,e) = [nElemNodes{1}(ex); nElemNodes{2}(ey)] - 1;
            if e == 1
                xIndex = (1 : 1 + (nElemNodes{1}(ex)-1));
                yIndex = (1 : 1 + (nElemNodes{2}(ey)-1));
            elseif ex == 1
                xIndex = (1 : 1 + (nElemNodes{1}(ex)-1));
                yIndex = ((ey-1)+1) + sum(nElemNodes{2}(1:ey-1)-2);
                yIndex = (yIndex(end) : (yIndex(end) + (nElemNodes{2}(ey)-1)));
            else
                xIndex = (xIndex(end) : (xIndex(end) + (nElemNodes{1}(ex)-1)));
            end
            nodeIDs = (yIndex-1)*nxNodes + xIndex';
            eCONN(1:numel(nodeIDs),e) = nodeIDs(:);
        end
    end
elseif nDimensions == 3
    nTotElem = prod(nElem);
    maxElemNodes = kron(kron(nElemNodes{1},nElemNodes{2}'),nElemNodes{3});
    maxElemNodes = max(maxElemNodes(:));
    eCONN = nan(maxElemNodes,nTotElem);
    
    nxNodes = (nElem(1)+1) + sum(nElemNodes{1}-2);
    nyNodes = (nElem(2)+1) + sum(nElemNodes{2}-2);
    nzNodes = (nElem(3)+1) + sum(nElemNodes{3}-2);
    
    NX = nElem(1);
    NY = nElem(2);
    NZ = nElem(3);
    for ez = 1:NZ
        for ey = 1:NY
            for ex = 1:NX
                e = (ez-1)*NY*NX + (ey-1)*NY + ex;
                eDegree(:,e) = [nElemNodes{1}(ex); nElemNodes{2}(ey); nElemNodes{3}(ez)] - 1;
                if ex == 1 && ey == 1 && ez == 1
                    xIndex = (1 : 1 + (nElemNodes{1}(ex)-1));
                    yIndex = (1 : 1 + (nElemNodes{2}(ey)-1));
                    zIndex = (1 : 1 + (nElemNodes{3}(ez)-1));
                elseif ex == 1 && ey == 1
                    xIndex = (1 : 1 + (nElemNodes{1}(ex)-1));
                    yIndex = (1 : 1 + (nElemNodes{2}(ey)-1));
                    zIndex = ((ez-1)+1) + sum(nElemNodes{1}(1:ex-1)-2) + sum(nElemNodes{2}(1:ey-1)-2);
                    zIndex = (zIndex(end) : (zIndex(end) + (nElemNodes{3}(ez)-1)));
                elseif ex == 1
                    xIndex = (1 : 1 + (nElemNodes{1}(ex)-1));
                    yIndex = ((ey-1)+1) + sum(nElemNodes{2}(1:ey-1)-2);
                    yIndex = (yIndex(end) : (yIndex(end) + (nElemNodes{2}(ey)-1)));
                else
                    xIndex = (xIndex(end) : (xIndex(end) + (nElemNodes{1}(ex)-1)));
                end
                zIndex = reshape(zIndex,[1,1,length(zIndex)]);
                nodeIDs = (zIndex-1)*nxNodes*nyNodes + (yIndex-1)*nxNodes + xIndex';
                eCONN(1:numel(nodeIDs),e) = nodeIDs(:);
            end
        end
    end
end
end