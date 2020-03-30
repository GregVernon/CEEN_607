clear

%% Load mesh
EXO = load("mesh/test_1x1.mat");

%% Initialize Mesh class
MESH = feMesh(EXO.nelems);

%% Initialize Elements
eID = 0;
nBlk = size(EXO.element_blocks,2);
for blk = 1:nBlk
    nBlkElems = size(EXO.element_blocks{4,blk},2);
    for e = 1:nBlk
        eID = eID + 1;
        MESH.Elements(eID).Degree = 1;
        MESH.Elements(eID).Dimension = 2;
        MESH.Elements(eID).NodeConnectivity = EXO.element_blocks{4,blk}(:,e);
        MESH.Configuration_Parametric = feElementConfiguration("Parametric");
        MESH.Configuration_Reference = feElementConfiguration("Reference");
    end
end
