include("feHeader.jl")

function feStructuralMechanics(inputFile)
    GEOM,PARAMS = importInput(inputFile);
    GEOM.Elements = buildElementQuadrature(GEOM.Elements, GEOM.Nodes);
    
    D̃ = computeMaterialStiffnessMatrix(1.,0.);

    GEOM.Elements = computeLocalElementStiffnessMatrices(GEOM.Elements,GEOM.Nodes,D̃);
    # K̃ = u -> assembleGlobalStiffnessMatrix(GEOM.Elements, GEOM.Nodes);
    
    # R̃ = u -> computeResidual(u, GEOM);
    # Ĩ = α -> InitialConditions(α, R̃, ũ, GEOM)
    # Ã = (R̃, K̃, u) -> applyBoundaryConditions(R̃, K̃, u, GEOM);
    # uₛ = newton_raphson(R̃, K̃, Ã, zeros(size(K̃(1),1)), 10, 10, 1e-5)
    uₛ = newton_raphson(GEOM, 10, 10, 1e-5)
    return uₛ
end