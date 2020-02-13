module feMesh

import ForwardDiff

include("feDatastruct.jl")
include("feQuadrature.jl")
include("feBasisFunctions.jl")
include("feEnumerations.jl")
import .feDatastruct
import .feQuadrature
import .feBasisFunctions
import .feEnumerations

export buildGlobalNodeCoordinateArray
export buildGlobalElementNodalConnectivityArray
export buildNodeGlobalDOFConnectivityArray
export buildConstrainedDOFList



function buildGlobalNodeCoordinateArray(G)
    # Get global information about the Genesis file
    num_dim       = G.dim["num_dim"]
    num_nodes     = G.dim["num_nodes"]
    # Load Nodal Positions
    if num_dim == 1
        coordx = G["coordx"].var[:]
        NodeCoords = coordx
    elseif num_dim == 2
        coordx = G["coordx"].var[:]
        coordy = G["coordy"].var[:]
        NodeCoords = [coordx coordy]
    elseif num_dim == 3
        coordx = G["coordx"].var[:]
        coordy = G["coordy"].var[:]
        coordz = G["coordz"].var[:]
        NodeCoords = [coordx coordy coordz]
    end
    
    return NodeCoords
end

function buildGlobalElementNodalConnectivityArray(ELEMS)
    num_elems = length(ELEMS)
    # Find max number of nodes in the elements
    max_loc_nodes = 0
    for e = 1:num_elems
        num_loc_nodes = ELEMS[e].NumNodes
        if num_loc_nodes > max_loc_nodes
            max_loc_nodes = ELEMS[e].NumNodes
        end
    end

    ElemConnect = fill(-1,max_loc_nodes,num_elems)
    for e = 1:num_elems
        num_loc_nodes = ELEMS[e].NumNodes
        ElemConnect[1:num_loc_nodes,e] = ELEMS[e].ChildNodes
    end

    return ElemConnect
end

function buildNodeGlobalDOFConnectivityArray(NODES)
    num_nodes = length(NODES)
    # Find max number of dofs in the nodes
    max_num_dofs = 0
    for n = 1:num_nodes
        nodeLocalDOFS = NODES[n].ChildDOFS
        num_local_dofs = length(nodeLocalDOFS)
        if num_local_dofs > max_num_dofs
            max_num_dofs = num_local_dofs
        end
    end

    # Populate the array
    NodeConnect = fill(-1,max_num_dofs,num_nodes)
    for n = 1:num_nodes
        nodeLocalDOFS = NODES[n].ChildDOFS
        num_local_dofs = length(nodeLocalDOFS)
        for ldof = 1:num_local_dofs
            NodeConnect[ldof,n] = nodeLocalDOFS[ldof]
        end
    end
    return NodeConnect
end

function buildConstrainedDOFList(NODES, NS)
    NodeConnect = buildNodeGlobalDOFConnectivityArray(NODES)
    num_gdof = max(NodeConnect...)
    isConstrainedGDOF = fill(false,num_gdof)
    num_nodesets = length(NS)
    for i = 1:num_nodesets
        constrained_ldof = NS[i].ConstrainedDOF
        num_ns_nodes = length(NS[i].ChildNodes)
        for n = 1:num_ns_nodes
            for ldof = 1:length(constrained_ldof)
                gdofID = NodeConnect[constrained_ldof[ldof],NS[i].ChildNodes[n]]
                isConstrainedGDOF[gdofID] = true
            end
        end
    end

    return isConstrainedGDOF
end

function buildElementQuadrature(ELEMS)
    num_elems = length(ELEMS)
    for e = 1:num_elems
        eDegree = ELEMS[e].Degree
        num_variates = ELEMS[e].NumVariates
        num_dims = ELEMS[e].Dimension
        ELEMS[e].Quadrature = feDatastruct.Quadrature()
        ELEMS[e].Quadrature.Type = "Gauss-Legendre"
        ELEMS[e].Quadrature.Component_Points = Array{Array{Float64},1}(undef,num_variates) #Array{Any,1}(undef,num_dims)
        ELEMS[e].Quadrature.Component_Weights = Array{Array{Float64},1}(undef,num_variates) #Array{Any,1}(undef,num_dims)
        nPts = zeros(Int64,num_variates)
        for dim = 1:num_variates
            nPts[dim] = Int(ceil((eDegree[dim] + 1)/2)) + 1
            ξ, W = feQuadrature.computeGaussQuadrature(nPts[dim])
            ELEMS[e].Quadrature.Component_Points[dim] = ξ
            ELEMS[e].Quadrature.Component_Weights[dim] = W
        end

        num_quad_points = prod(nPts)
        ELEMS[e].Quadrature.Points = fill(zeros(num_dims),num_quad_points)
        ELEMS[e].Quadrature.Weights = zeros(num_quad_points)
        q = 0
        if num_variates == 1
            for q1 = 1:nPts[1]
                q+=1
                ELEMS[e].Quadrature.Points[q] = ELEMS[e].Quadrature.Component_Points[1][q1]
                ELEMS[e].Quadrature.Weights[q] = ELEMS[e].Quadrature.Component_Weights[1][q1]
            end
        elseif num_variates == 2
            for q1 = 1:nPts[1]
                for q2 = 1:nPts[2]
                    q+=1
                    ELEMS[e].Quadrature.Points[q] = [ELEMS[e].Quadrature.Component_Points[1][q1], ELEMS[e].Quadrature.Component_Points[2][q2]]
                    ELEMS[e].Quadrature.Weights[q] = ELEMS[e].Quadrature.Component_Weights[1][q1] * ELEMS[e].Quadrature.Component_Weights[2][q2]
                end
            end
        elseif num_variates == 3
            for q1 = 1:nPts[1]
                for q2 = 1:nPts[2]
                    for q3 = 1:nPts[3]
                        q+=1
                        ELEMS[e].Quadrature.Points[q] = [ELEMS[e].Quadrature.Component_Points[1][q1], ELEMS[e].Quadrature.Component_Points[2][q2], ELEMS[e].Quadrature.Component_Points[3][q3]]
                        ELEMS[e].Quadrature.Weights[q] = ELEMS[e].Quadrature.Component_Weights[1][q1] * ELEMS[e].Quadrature.Component_Weights[2][q2] * ELEMS[e].Quadrature.Component_Weights[3][q3]
                    end
                end
            end
        end
    end
    return ELEMS
end

function buildElementBasisFunctions(ELEMS)
    num_elems = length(ELEMS)
    for e = 1:num_elems
        elem_degree = ELEMS[e].Degree
        num_dims = ELEMS[e].Dimension
        num_loc_nodes = ELEMS[e].NumNodes
        bfun = Array{Any,1}(undef,num_loc_nodes)
        for n = 1:num_loc_nodes
            bfun[n] = feBasisFunctions.constructBasis(Val(feBasisFunctions.feEnumerations.lagrange),elem_degree,n)
        end
        ELEMS[e].Basis = bfun
    end
    return ELEMS
end

function buildElementBasisGradientFunctions(ELEMS)
    num_elems = length(ELEMS)
    for e = 1:num_elems
        currElement = ELEMS[e]
        elem_degree = currElement.Degree
        elem_dim = currElement.Dimension
        elemBasis = currElement.Basis
        num_dims = length(elem_degree)
        num_loc_nodes = currElement.NumNodes
        gfun = Array{Any,1}(undef,num_loc_nodes)
        for n = 1:num_loc_nodes
            gfun[n] = ξ->ForwardDiff.gradient(elemBasis[n],ξ)
        end
        ELEMS[e].∂Basis = gfun
    end
    return ELEMS
end

function cacheQuadratureBasisEvaluations(ELEMS)
    num_elems = length(ELEMS)
    for e = 1:num_elems
        num_quad_points = length(ELEMS[e].Quadrature.Points)
        num_loc_nodes = ELEMS[e].NumNodes
        num_variates = ELEMS[e].NumVariates
        num_dim = ELEMS[e].Dimension
        ELEMS[e].Quadrature.Basis_Evaluation = zeros(num_quad_points)
        ELEMS[e].Quadrature.∂Basis_Evaluation = fill(zeros(num_dim),num_loc_nodes)
        for q = 1:num_quad_points
            ξ = zeros(num_dim)
            ξ[1:num_variates] .= ELEMS[e].Quadrature.Points[q]
            for n = 1:num_loc_nodes
                ELEMS[e].Quadrature.Basis_Evaluation[q] += ELEMS[e].Basis[n](ξ)
                ELEMS[e].Quadrature.∂Basis_Evaluation[q] +=  ELEMS[e].∂Basis[n](ξ)
            end
        end
    end
    return ELEMS
end

function cacheQuadratureJacobianEvaluations(ELEMS,NODES)
    num_elems = length(ELEMS)
    for e = 1:num_elems
        num_quad_points = length(ELEMS[e].Quadrature.Points)
        num_local_nodes = ELEMS[e].NumNodes
        num_variates = ELEMS[e].NumVariates
        num_dim = ELEMS[e].Dimension
        ELEMS[e].Quadrature.JacobianMatrix_Evaluation = Array{Array{Float64,2},1}(undef,num_quad_points)
        for q = 1:num_local_nodes
            ELEMS[e].Quadrature.JacobianMatrix_Evaluation[q] = zeros(Float64,num_local_nodes,num_dim)
        end

        for n = 1:num_local_nodes
            x = NODES[ELEMS[e].ChildNodes[n]].Coordinates
            for q = 1:num_quad_points
                ξ = zeros(num_dim)
                ξ[1:num_variates] .= ELEMS[e].Quadrature.Points[q]
                J = computeJacobian(ELEMS[e],ξ,x)
                ELEMS[e].Quadrature.JacobianMatrix_Evaluation[q][1:size(J,1),1:size(J,2)] .+= J
            end
        end
    end
    return ELEMS
end

function Parametric_2_Cartesian(Element,ξ)
    num_loc_nodes = Element.NumNodes
    num_dim = length(ξ)
    x = zeros(num_dim)
    for n = 1:num_loc_nodes
        x += Element.Basis[n](ξ)
    end
    return x
end

function computeJacobian(Element,ξ,x)
    num_loc_nodes = Element.NumNodes
    num_variates = Element.NumVariates
    num_dim = Element.Dimension
    J = zeros(Float64,num_loc_nodes,num_dim)
    for n = 1:num_loc_nodes
        J[n,:] = Element.∂Basis[n](ξ) .* x
    end
    return J
end

function getNumVariates(elem_type::String)
    ETYPE_1D = ["BAR2"]
    ETYPE_2D = ["TRI3","QUAD4"]
    ETYPE_3D = ["TETRA4","PYRAMID5","WEDGE6","HEX8"]
    if      elem_type in ETYPE_1D
        NumVar = 1
    elseif  elem_type in ETYPE_2D
        NumVar = 2
    elseif  elem_type in ETYPE_3D
        NumVar = 3
    end
    return NumVar
end

function assignNodeDOFS(NODES, num_dim)
    dofID = 0
    for n = 1:length(NODES)
        NODES[n].ChildDOFS = zeros(Int64,num_dim)
        for d = 1:num_dim
            dofID+=1
            NODES[n].ChildDOFS[d] = dofID
        end
    end
    return NODES
end

function setParentElement(ELEMS,NODES)
    num_nodes = length(NODES)
    for n = 1:num_nodes
        NODES[n].ParentElements = Int64[]
    end

    num_elem = length(ELEMS)
    for e = 1:num_elem
        num_nodes = length(ELEMS[e].ChildNodes)
        for n = 1:num_nodes
            gnID = ELEMS[e].ChildNodes[n]
            if (e in NODES[gnID].ParentElements) == false
                append!(NODES[gnID].ParentElements,e)
            end
        end
    end
    return NODES
end

function setElementNodeTypes(ELEMS,NODES)
    num_elem = length(ELEMS)
    for e = 1:num_elem
        num_nodes = length(ELEMS[e].ChildNodes)
        ELEMS[e].BoundaryNodes = zeros(Int64,num_nodes)
        ELEMS[e].CornerNodes = zeros(Int64,num_nodes)
        ELEMS[e].FaceNodes = zeros(Int64,num_nodes)
        ELEMS[e].InternalNodes = zeros(Int64,num_nodes)
        for n = 1:num_nodes
            gnID = ELEMS[e].ChildNodes[n]
            if NODES[gnID].isElementBoundaryNode
                ELEMS[e].BoundaryNodes[n] = gnID
            end
            if NODES[gnID].isElementCornerNode
                ELEMS[e].CornerNodes[n] = gnID
            end
            if NODES[gnID].isElementFaceNode
                ELEMS[e].FaceNodes[n] = gnID
            end
            if NODES[gnID].isElementInternalNode
                ELEMS[e].InternalNodes[n] = gnID
            end
        end
    end
    
    return ELEMS
end

end