%% 2D
clear
xMin = 0;
xMax = 1;
yMin = 0;
yMax = 1;

nElem{1} = 4;
eDegree{1} = [1 2 3 4];
eDomains{1} = linspace(xMin,xMax,nElem{1}+1);
nElemNodes{1} = eDegree{1}+1;

nElem{2} = 2;
eDegree{2} = [4 3];
eDomains{2} = linspace(yMin,yMax,nElem{2}+1);
nElemNodes{2} = eDegree{2}+1;

% Mesh the 2D nodes
[x,y] = meshNodes(eDomains,nElemNodes);
figure
scatter(x,y,'filled')
hold on
nodeID = 0;
for ii = 1:length(x)
    nodeID = nodeID+1;
    text(x(ii),y(ii),num2str(nodeID),'FontWeight','Bold','FontSize',16)
end

% Build element connectivity
[eCONN, eDEGREE] = buildConnectivity(nElemNodes);

GNode2GDOF = buildGlobalNodeDOFS(length(x),2);

globalBoundaryNodes = getGlobalBoundaryNodeIDs(eCONN,eDEGREE);

scatter(x(globalBoundaryNodes),y(globalBoundaryNodes),'filled')