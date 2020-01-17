function varargout = meshNodes(eDomains,nElemNodes)
    nDimensions = length(eDomains);
    nDimensions = 1;

if nDimensions == 1
    nElem = length(nElemNodes{1});
    nNodes = (nElem+1) + sum(nElemNodes{1}-2);
    x = zeros(nNodes,1);
    for e = 1:nElem
        if e == 1
            nodeIDs = 1:nElemNodes{1}(e);
        else
            nodeIDs = nodeIDs(end):(nodeIDs(end)+(nElemNodes{1}(e)-1));
        end

        x(nodeIDs) = linspace(eDomains{1}(e),eDomains{1}(e+1),nElemNodes{1}(e));
    end
    varargout{1} = x;
elseif nDimensions == 2
    x = meshNodes(eDomains{1},nElemNodes{1});
    y = meshNodes(eDomains{2},nElemNodes{2});

    XX = zeros(length(x)*length(y),1);
    YY = zeros(length(x)*length(y),1);
    for jj = 1:length(y)
        for ii = 1:length(x)
            nodeID = (jj-1)*length(x) + ii;
            XX(nodeID) = x(ii);
            YY(nodeID) = y(jj);
        end
    end
    varargout{1} = XX;
    varargout{2} = YY;
    
elseif nDimensions == 3
    x = meshNodes(eDomains{1},nElemNodes{1});
    y = meshNodes(eDomains{2},nElemNodes{2});
    z = meshNodes(eDomains{3},nElemNodes{3});

    XXX = zeros(length(x)*length(y)*length(z),1);
    YYY = zeros(length(x)*length(y)*length(z),1);
    ZZZ = zeros(length(x)*length(y)*length(z),1);
    for kk = 1:length(z)
        for jj = 1:length(y)
            for ii = 1:length(x)
                nodeID = (kk-1)*length(x)*length(y) + (jj-1)*length(x) + ii;
                XXX(nodeID) = x(ii);
                YYY(nodeID) = y(jj);
                ZZZ(nodeID) = z(kk);
            end
        end
    end
    varargout{1} = XXX;
    varargout{2} = YYY;
    varargout{3} = ZZZ;
end
end
