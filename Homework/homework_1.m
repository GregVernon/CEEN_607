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
% Build element-node connectivity
eCONN = buildConnectivity_1D(nElemNodes);

plot(x,zeros(size(x)),'-o')
disp(eCONN)

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

nElem{2} = 4;
eDegree{2} = [4 3 2 1];
eDomains{2} = linspace(yMin,yMax,nElem{2}+1);
nElemNodes{2} = eDegree{2}+1;

% Mesh the 2D nodes
[x,y] = meshNodes_2D(eDomains,nElemNodes);
% figure
% scatter(x,y,'filled')
% hold on
% nodeID = 0;
% for ii = 1:length(x)
%     nodeID = nodeID+1;
%     text(x(ii),y(ii),num2str(nodeID),'FontWeight','Bold','FontSize',16)
% end

% Build element connectivity
eCONN = buildConnectivity_2D(nElemNodes);

% figure
% hold on
% cmap = lines(size(eCONN,2));
% for e = 1:size(eCONN,2)
%     elemNodes = ~isnan(eCONN(:,e));
%     px = [x(eCONN(elemNodes,e)); x(eCONN(1,e))];
%     py = [y(eCONN(elemNodes,e)); y(eCONN(1,e))];
%     patch(px,py,cmap(e,:))
% end
% scatter(x,y,'filled')

%% 3D
clear
xMin = 0;
xMax = 1;
yMin = 0;
yMax = 1;
zMin = 0;
zMax = 1;

nElem{1} = 4;
eDegree{1} = [1 2 3 4];
eDomains{1} = linspace(xMin,xMax,nElem{1}+1);
nElemNodes{1} = eDegree{1}+1;

nElem{2} = 4;
eDegree{2} = [4 3 2 1];
eDomains{2} = linspace(yMin,yMax,nElem{2}+1);
nElemNodes{2} = eDegree{2}+1;

nElem{3} = 4;
eDegree{3} = [1 2 2 1];
eDomains{3} = linspace(zMin,zMax,nElem{3}+1);
nElemNodes{3} = eDegree{3}+1;

% Mesh the 3D nodes
[x,y,z] = meshNodes_3D(eDomains,nElemNodes);

% Build element connectivity
eCONN = buildConnectivity_3D(nElemNodes);

%% Mesh Generation Functions
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

function [XX,YY] = meshNodes_2D(eDomains,nElemNodes)
x = meshNodes_1D(eDomains{1},nElemNodes{1});
y = meshNodes_1D(eDomains{2},nElemNodes{2});

XX = zeros(length(x)*length(y),1);
YY = zeros(length(x)*length(y),1);
for jj = 1:length(y)
    for ii = 1:length(x)
        nodeID = (jj-1)*length(x) + ii;
        XX(nodeID) = x(ii);
        YY(nodeID) = y(jj);
    end
end
end

function [XXX,YYY,ZZZ] = meshNodes_3D(eDomains,nElemNodes)
x = meshNodes_1D(eDomains{1},nElemNodes{1});
y = meshNodes_1D(eDomains{2},nElemNodes{2});
z = meshNodes_1D(eDomains{3},nElemNodes{3});

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
end


%% Element Connecitivty Functions
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

function eCONN = buildConnectivity_2D(nElemNodes)
nElem = cellfun(@length,nElemNodes);
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
end

function eCONN = buildConnectivity_3D(nElemNodes)
nElem = cellfun(@length,nElemNodes);
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
            e = (ey-1)*NY + ex;
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