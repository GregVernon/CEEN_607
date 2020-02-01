%% 1D
clear
xMin = 0;
xMax = 1;
yMin = 0;
yMax = 1;

nElem{1} = 4;
eDegree{1} = [1 2 3 4];
eDomains{1} = linspace(xMin,xMax,nElem{1}+1);
nElemNodes{1} = eDegree{1}+1;

% Mesh the 1D nodes
x = meshNodes(eDomains,nElemNodes);
figure
scatter(x,zeros(size(x)),'filled')
hold on
nodeID = 0;
for ii = 1:length(x)
    nodeID = nodeID+1;
    text(x(ii),1,num2str(nodeID),'FontWeight','Bold','FontSize',16)
end
ylim([-10 10])

% Build element connectivity
[eCONN, eDEGREE] = buildConnectivity(nElemNodes);

GNode2GDOF = buildGlobalNodeDOFS(length(x),2);

globalBoundaryNodes = getGlobalBoundaryNodeIDs(eCONN,eDEGREE);

scatter(x(globalBoundaryNodes),zeros(size(globalBoundaryNodes)),'filled')