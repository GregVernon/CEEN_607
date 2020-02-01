%% 3D
clear
xMin = 0;
xMax = 1;
yMin = 0;
yMax = 1;
zMin = 0;
zMax = 1;

nElem{1} = 4;
eDegree{1} = [4 3 2 1];
eDomains{1} = linspace(xMin,xMax,nElem{1}+1);
nElemNodes{1} = eDegree{1}+1;

nElem{2} = 4;
eDegree{2} = [1 2 3 4];
eDomains{2} = linspace(yMin,yMax,nElem{2}+1);
nElemNodes{2} = eDegree{2}+1;

nElem{3} = 4;
eDegree{3} = [1 2 2 1];
eDomains{3} = linspace(zMin,zMax,nElem{3}+1);
nElemNodes{3} = eDegree{3}+1;

% Mesh the 3D nodes
[x,y,z] = meshNodes(eDomains,nElemNodes);

% Build element connectivity
[eCONN, eDEGREE] = buildConnectivity(nElemNodes);

figure
scatter3(x,y,z,'filled')
hold on
nodeID = 0;
for ii = 1:length(x)
    nodeID = nodeID+1;
    text(x(ii),y(ii),z(ii),num2str(nodeID),'FontWeight','Bold','FontSize',8)
end
axis vis3d
view(60,30)

GNode2GDOF = buildGlobalNodeDOFS(length(x),2);

globalBoundaryNodes = getGlobalBoundaryNodeIDs(eCONN,eDEGREE);

scatter3(x(globalBoundaryNodes),y(globalBoundaryNodes),z(globalBoundaryNodes),'filled')