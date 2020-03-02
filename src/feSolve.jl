include("feHeader.jl")

function feStructuralMechanics(inputFile)
    GEOM,PARAMS = importInput(inputFile)
    GEOM.Elements = buildElementQuadrature(GEOM.Elements)
    
    D̃ = computeMaterialStiffnessMatrix(1000,0.1)

    GEOM.Elements = computeLocalElementStiffnessMatrices(GEOM.Elements,GEOM.Nodes,D̃)
    K̃ = u -> assembleGlobalStiffnessMatrix(GEOM.Elements, GEOM.Nodes)
    
    R̃ = u -> computeResidual(u, GEOM)
    uₛ = newton_raphson(R̃, K̃, zeros(size(K̃(1),1)), 10, 10, 1e-5)
end