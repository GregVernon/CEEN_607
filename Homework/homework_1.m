%% 1D
clear
xMin = 0;
xMax = 1;

nElem = 4;
eDegree = [1 2 3 4];
eDomains = linspace(xMin,xMax,nElem+1);
nElemNodes = eDegree+1;

% Mesh the 1D nodes
x = meshNodes_1D(eDomains,nElemNodes);
eCONN = buildConnectivity_1D(nElemNodes);

plot(x,zeros(size(x)),'-o')
disp(eCONN)

%% Functions
function x = meshNodes_1D(eDomains,nElemNodes)
nElem = length(nElemNodes);
nNodes = (nElem+1) + sum(nElemNodes-2);
x = zeros(nNodes,1);
for e = 1:nElem
    if e == 1
        nodeIDs = 1:nElemNodes(e);
    else
        nodeIDs = nodeIDs(end):(nodeIDs(end)+(nElemNodes(e)-1));
    end
    
    x(nodeIDs) = linspace(eDomains(e),eDomains(e+1),nElemNodes(e));
end
end

function eCONN = buildConnectivity_1D(nElemNodes)
nElem = length(nElemNodes);
eCONN = nan(max(nElemNodes),nElem);
for e = 1:nElem
    if e == 1
        nodeIDs = 1:nElemNodes(e);
    else
        nodeIDs = nodeIDs(end):(nodeIDs(end)+(nElemNodes(e)-1));
    end
    eCONN(1:nElemNodes(e),e) = nodeIDs;
end
end