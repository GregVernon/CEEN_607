using NamedDims
import LinearAlgebra
import ForwardDiff

function GaussQuadrature_1D(fun, nPts)
    ξ,W = GaussQuadratureRule_1D(nPts)
    nQP = nPts^1

    ∫fun = 0.
    for n = 1:nQP
        ∫fun += fun(ξ[n]) * W[n]
    end
    return ∫fun
end

function GaussQuadrature_2D(fun,nPts)
    ξ,W = GaussQuadratureRule_2D(nPts)
    nQP = nPts^2

    sizeFun = size(fun(ξ[1,:]))
    if sizeFun == ()
        ∫fun = 0.
    else
        ∫fun = zeros(Float64, sizeFun)
    end

    for n = 1:nQP
        ∫fun += fun(ξ[n,:]) * W[n]
    end
    return ∫fun
end

function GaussQuadrature_3D(fun,nPts)
    ξ,W = GaussQuadratureRule_3D(nPts)
    nQP = nPts^3

    ∫fun = 0.
    for n = 1:nQP
        ∫fun += fun(ξ[n,:]) * W[n]
    end
    return ∫fun

end

function GaussQuadratureRule_1D(nPts)
    if nPts == 1
        ξ = NamedDimsArray{(:local_qp_id,)}([0.])
        W = NamedDimsArray{(:local_qp_id,)}([2.])
    elseif nPts == 2
        ξ = NamedDimsArray{(:local_qp_id,)}([-1/sqrt(3), 1/sqrt(3)])
        W = NamedDimsArray{(:local_qp_id,)}([1., 1.,])
    elseif nPts == 3
        ξ = NamedDimsArray{(:local_qp_id,)}([-sqrt(3)/sqrt(5), 0., sqrt(3)/sqrt(5)])
        W = NamedDimsArray{(:local_qp_id,)}([5/9, 8/9, 5/9])
    end
    return ξ,W
end

function GaussQuadratureRule_2D(nPts)
    W = Array{Array{Float64},1}(undef,2)
    ξ,W[1] = GaussQuadratureRule_1D(nPts)
    η,W[2] = GaussQuadratureRule_1D(nPts)

    nQP = nPts^2
    QP = NamedDimsArray{(:local_qp_id, :ℝᴺ)}(zeros(Float64,nQP,2))
    QW = NamedDimsArray{(:local_qp_id,)}(zeros(Float64,nQP))
    n = 0
    for j = 1:nPts
        for i = 1:nPts
            n+=1
            QP[n,:] = [ξ[i], η[j]]
            QW[n] = W[1][i] * W[2][j]
        end
    end
    return QP,QW
end

function GaussQuadratureRule_3D(nPts)
    W = Array{Array{Float64},1}(undef,3)
    ξ,W[1] = GaussQuadratureRule_1D(nPts)
    η,W[2] = GaussQuadratureRule_1D(nPts)
    ζ,W[3] = GaussQuadratureRule_1D(nPts)

    nQP = nPts^3
    QP = NamedDimsArray{(:local_qp_id, :ℝᴺ)}(zeros(Float64,nQP,3))
    QW = NamedDimsArray{(:local_qp_id,)}(zeros(Float64,nQP))
    n = 0
    for k = 1:nPts
        for j = 1:nPts
            for i = 1:nPts
                n+=1
                QP[n,:] = [ξ[i], η[j], ζ[k]]
                QW[n] = W[1][i] * W[2][j] * W[3][k]
            end
        end
    end
    return QP,QW
end

function LagrangeBasis_1D(deg, ξ)
    ξ = ξ[1]
    L = NamedDimsArray{(:local_node_id,)}(Array{Any,1}(undef,deg+1))
    if deg == 1
        L = [(1-ξ)/2, (1+ξ)/2]
        return L
    elseif deg == 2
        L = [(ξ^2 - ξ)/2, 1-ξ^2, (ξ^2 + ξ)/2]
        return L
    end
end

function LagrangeBasis_2D(deg,ξ)
    L = NamedDimsArray{(:local_node_id,)}(Array{Any,1}(undef,(deg+1)^2))
    L1 = LagrangeBasis_1D(deg,ξ[1])
    L2 = LagrangeBasis_1D(deg,ξ[2])
    n = 0
    for j = 1:deg+1
        for i = 1:deg+1
            n += 1
            L[n] = L1[i] * L2[j]
        end
    end
    return L
end

function LagrangeBasis_3D(deg,ξ)
    L = NamedDimsArray{(:local_node_id,)}(Array{Any,1}(undef,(deg+1)^3))
    L1 = LagrangeBasis_1D(deg,ξ[1])
    L2 = LagrangeBasis_1D(deg,ξ[2])
    L3 = LagrangeBasis_1D(deg,ξ[3])
    n = 0
    for k = 1:deg+1
        for j = 1:deg+1
            for i = 1:deg+1
                n += 1
                L[n] = L1[i] * L2[j] * L3[k]
            end
        end
    end
    return L
end

function ∇LagrangeBasis_1D(deg,ξ)
    f(ξ) = LagrangeBasis_1D(deg,ξ)
    ∇f = ForwardDiff.jacobian(ξ->f(ξ),ξ)
    return ∇f
end

function ∇LagrangeBasis_2D(deg,ξ)
    f(ξ) = LagrangeBasis_2D(deg,ξ)
    ∇f = ForwardDiff.jacobian(ξ->f(ξ),ξ)
    return ∇f
end

function ∇LagrangeBasis_3D(deg,ξ)
    f(ξ) = LagrangeBasis_3D(deg,ξ)
    ∇f = ForwardDiff.jacobian(ξ->f(ξ),ξ)
    return ∇f
end

function computeGeometricMapping(Nₐ, xₐ, ξ)
    x = NamedDimsArray{(:ℝᴺ,)}(zeros(Float64,size(ξ)))
    num_nodes = size(xₐ, :local_node_id)
    for n = 1:num_nodes
        x += Nₐ(ξ)[n] .* xₐ[n,:]
    end
    return x
end

function compute∇GeometricMapping(∇Nₐ, xₐ, ξ)
    num_nodes = size(xₐ,:local_node_id)
    num_par_dim = length(ξ)
    num_cart_dim = size(xₐ,:ℝᴺ)
    Jᵢⱼ = NamedDimsArray{(:ℙᴺ,:ℝᴺ,)}(zeros(Float64,length(ξ),num_cart_dim))

    for n = 1:num_nodes
        for j = 1:num_par_dim
            for i = 1:num_cart_dim
                Jᵢⱼ[j,i] += ∇Nₐ(ξ)[n,j] * xₐ[n,i] 
            end
        end
    end
    return transpose(Jᵢⱼ)
end

function compute∇ₓNₐ(∇Nₐ, Jᵢⱼ)
    ∇ₓNₐ = ∇Nₐ * inv(Jᵢⱼ)
    return ∇ₓNₐ
end

function StrainDisplacement_2D(∇Nₐ)
    B = [∇Nₐ[1]  0.0;
         0.0    ∇Nₐ[2];
         ∇Nₐ[2] ∇Nₐ[1]]
end

function computeVirtualStrain_2D(∇Nₐ, cₐ, ξ)
    num_nodes = length(cₐ)
    εᵥ = zeros(Float64,3)
    for a = 1:num_nodes
        εᵥ += StrainDisplacement_2D(∇Nₐ[a,:]) * cₐ
    end
    return εᵥ
end

function computeIntegralScaling_2D(Jᵢⱼ, sideID)
    if sideID == 0
        scaleFactor = LinearAlgebra.det(Jᵢⱼ)  
    elseif sideID == 1 || sideID == 2
        scaleFactor = LinearAlgebra.norm(Jᵢⱼ[:,2])
    elseif sideID == 3 || sideID == 4
        scaleFactor = LinearAlgebra.norm(Jᵢⱼ[:,1])
    end
    return scaleFactor
end

function computeBoundaryNormals(Jᵢⱼ, sideID)
    if sideID == 0
        ñ = undef
    elseif sideID == 1
        ñ = 1/LinearAlgebra.norm(LinearAlgebra.cross([0,0,1], [Jᵢⱼ[:,2]..., 0.0])) * (LinearAlgebra.cross([0,0,1], [Jᵢⱼ[:,2]..., 0.0]))
    elseif sideID == 2
        ñ = 1/LinearAlgebra.norm(LinearAlgebra.cross([Jᵢⱼ[:,2]..., 0.0], [0,0,1])) * (LinearAlgebra.cross([Jᵢⱼ[:,2]..., 0.0], [0,0,1]))
    elseif sideID == 3
        ñ = 1/LinearAlgebra.norm(LinearAlgebra.cross([Jᵢⱼ[:,1]..., 0.0], [0,0,1])) * (LinearAlgebra.cross([Jᵢⱼ[:,1]..., 0.0], [0,0,1]))
    elseif sideID == 4
        ñ = 1/LinearAlgebra.norm(LinearAlgebra.cross([0,0,1], [Jᵢⱼ[:,1]..., 0.0])) * (LinearAlgebra.cross([0,0,1], [Jᵢⱼ[:,1]..., 0.0]))
    end
    return ñ[1:2]
end

function computeMaterialStiffnessMatrix(E, ν)
    D̃ = [1-ν   ν    0;
          ν   1-ν   0;
          0    0   1-2ν]
    
    D̃ *= E/((1+ν)*(1-2ν))
    return D̃
end

function computeExternalForce(GEOM)
    GEOM.Elements = computeLocalExternalForceVector(GEOM.Elements, GEOM.Nodes, GEOM.ElementSets, GEOM.SurfaceSets, GEOM.NodeSets)
    F_ext = assembleGlobalExternalForceVector(GEOM.Elements, GEOM.Nodes)
end

function computeInternalForce(Δu, GEOM)
    GEOM.Elements = computeLocalInternalForceVector(Δu, GEOM.Elements, GEOM.Nodes)
    F_int = assembleGlobalInternalForceVector(GEOM.Elements, GEOM.Nodes)
end

function computeResidual(Δu, GEOM)
    F_External = computeExternalForce(GEOM)
    F_Internal = computeInternalForce(Δu, GEOM)
    Residual = F_External - F_Internal
    return Residual
end

function computeLocalInternalForceVector(Δu, ELEMS, NODES)
    num_elem = length(ELEMS)
    num_dof_per_node = length(NODES[1].ChildDOFS)

    D̃ = computeMaterialStiffnessMatrix(1000, 0.1)

    # Initialize each element's local internal force vector to zero-vector
    for e = 1:num_elem
        num_loc_nodes = length(ELEMS[e].ChildNodes)
        ELEMS[e].InternalForceVector = NamedDimsArray{(:local_dof_id,)}(zeros(Float64, num_loc_nodes*num_dof_per_node))
    end

    for e = 1:num_elem
        nPts = 2
        eDegree = ELEMS[e].Degree
        num_loc_nodes = length(ELEMS[e].ChildNodes)
        side_id = 1
        for qp_id = 1:size(ELEMS[e].Quadrature[1].QuadraturePoints,:local_qp_id)
            QP = ELEMS[e].Quadrature[side_id].QuadraturePoints[qp_id]
            ξ,W = GaussQuadratureRule_2D(nPts) # QP.Coordinates
            for n = 1:num_loc_nodes
                local_dofs = (n-1) * num_dof_per_node .+ collect(1:num_dof_per_node)
                global_dofs = NODES[ELEMS[e].ChildNodes[n]].ChildDOFS
                ∇Nₐ = ∇LagrangeBasis_2D(eDegree[1],QP.Coordinates)
                B̃ = StrainDisplacement_2D(∇Nₐ[n,:])
                ε = computeVirtualStrain_2D(∇Nₐ, Δu[global_dofs], QP.Coordinates)
                σ =  D̃ * ε
                ELEMS[e].InternalForceVector[local_dofs] += transpose(B̃) * σ * LinearAlgebra.det(QP.Jᵢⱼ) * QP.Weights
            end
        end
    end
    
    return ELEMS
end

function assembleGlobalInternalForceVector(ELEMS, NODES)
    num_elem = length(ELEMS)
    num_nodes = length(NODES)
    num_dof_per_node = length(NODES[1].ChildDOFS)
    num_global_dofs = num_nodes * num_dof_per_node
    GlobalInternalForceVector = NamedDimsArray{(:global_dof_id,)}(zeros(Float64, num_global_dofs))
    for e = 1:num_elem
        num_loc_nodes = length(ELEMS[e].ChildNodes)
        for n = 1:num_loc_nodes
            local_dofs = (n -1) * num_dof_per_node .+ collect(1:num_dof_per_node)
            global_dofs = NODES[ELEMS[e].ChildNodes[n]].ChildDOFS
            GlobalInternalForceVector[global_dofs] += collect(ELEMS[e].InternalForceVector[local_dofs])
        end
    end
    return GlobalInternalForceVector
end

function computeLocalExternalForceVector(ELEMS, NODES, ELEMENTSETS, SURFACESETS, NODESETS)
    num_elem = length(ELEMS)
    num_dof_per_node = length(NODES[1].ChildDOFS)
    num_nodes = length(NODES)
    
    # Initialize each element's local external force vector to zero-vector
    for e = 1:num_elem
        num_loc_sides = size(ELEMS[e].SideNodes, :local_side_id)
        ELEMS[e].ExternalForceVector = NamedDimsArray{(:local_dof_id,)}(zeros(Float64, num_nodes * num_dof_per_node))
    end

    # Evaluate forces on element sets
    for es_id = 1:length(ELEMENTSETS)
        if isempty(ELEMENTSETS[es_id].LC_Type) == false
            for lc_id = 1:length(ELEMENTSETS[es_id].LC_Type)
                load_type = ELEMENTSETS[es_id].LC_Type[lc_id]
                if Int(load_type) == Int(feEnumerate.body)
                    f̃ = ELEMENTSETS[es_id].LC_Magnitude[lc_id] * ELEMENTSETS[es_id].LC_Direction / LinearAlgebra.norm(ELEMENTSETS[es_id].LC_Direction)
                    for loc_elem_id = 1:length(ELEMENTSETS[es_id].ChildElements)
                        global_elem_id = ELEMENTSETS[es_id].ChildElements[loc_elem_id]
                        local_side_id = 1
                        for qp_id = 1:length(ELEMS[global_elem_id].Quadrature[local_side_id].QuadraturePoints)
                            QP = ELEMS[global_elem_id].Quadrature[local_side_id].QuadraturePoints[qp_id]
                            for loc_node_id = 1:length(ELEMS[global_elem_id].ChildNodes)
                                loc_dof_id = (loc_node_id-1) * num_dof_per_node .+ collect(1:num_dof_per_node)
                                ELEMS[global_elem_id].ExternalForceVector[loc_dof_id] +=  f̃ * QP.Nₐ[loc_node_id] * LinearAlgebra.det(QP.Jᵢⱼ) * QP.Weights
                            end
                        end
                    end
                end
            end
        end
    end

    # Evaluate forces on surface sets
    for ss_id = 1:length(SURFACESETS)
        if isempty(SURFACESETS[ss_id].LC_Type) == false
            for lc_id = 1:length(SURFACESETS[ss_id].LC_Type)
                load_type = SURFACESETS[ss_id].LC_Type[lc_id]
                if Int(load_type) == Int(feEnumerate.pressure)
                    for loc_elem_id = 1:length(SURFACESETS[ss_id].ChildElements)
                        global_elem_id = SURFACESETS[ss_id].ChildElements[loc_elem_id]
                        side_id = SURFACESETS[ss_id].ChildElements_LocalFace[loc_elem_id]
                        for qp_id = 1:length(ELEMS[global_elem_id].Quadrature[side_id].QuadraturePoints)
                            QP = ELEMS[global_elem_id].Quadrature[side_id].QuadraturePoints[qp_id]
                            f̃ = SURFACESETS[ss_id].LC_Magnitude[lc_id] * QP.ñ / LinearAlgebra.norm(QP.ñ)
                            for loc_node_id = 1:length(ELEMS[global_elem_id].ChildNodes)
                                loc_dof_id = (loc_node_id-1) * num_dof_per_node .+ collect(1:num_dof_per_node)
                                ELEMS[global_elem_id].ExternalForceVector[loc_dof_id] +=  f̃ * QP.Nₐ[loc_node_id] * LinearAlgebra.det(QP.Jᵢⱼ) * QP.Weights
                            end
                        end
                    end
                end
            end
        end
    end

    # Evaluate forces on node sets
    for ns_id = 1:length(NODESETS)
        if isdefined(NODESETS[ns_id],:LC_Type) && isempty(NODESETS[ns_id].LC_Type) == false
            for lc_id = 1:length(NODESETS[ns_id].LC_Type)
                load_type = NODESETS[ns_id].LC_Type[lc_id]
                if Int(load_type) == Int(feEnumerate.force)
                    f̃ = NODESETS[ns_id].LC_
                    for ns_node_id = 1:length(NODESETS[ns_id].ChildNodes)
                        global_node_id = NODESETS[ns_id].ChildNodes
                        node_parent_elems = NODES[global_node_id].ParentElements
                        for e = 1:length(node_parent_elems)
                            loc_node_id = findall(ELEMS == global_node_id)
                            loc_dof_id = (loc_node_id -1) * num_dof_per_node .+ collect(1:num_dof_per_node)
                            ELEMS[node_parent_elems[e]].ExternalForceVector[loc_dof_id] += f̃
                        end
                    end
                end
            end
        end
    end

    return ELEMS
end

function assembleGlobalExternalForceVector(ELEMS, NODES)
    num_elem = length(ELEMS)
    num_nodes = length(NODES)
    num_dof_per_node = length(NODES[1].ChildDOFS)
    num_global_dofs = num_nodes * num_dof_per_node
    GlobalExternalForceVector = NamedDimsArray{(:global_dof_id,)}(zeros(Float64, num_global_dofs))
    for e = 1:num_elem
        num_loc_nodes = length(ELEMS[e].ChildNodes)
        for n = 1:num_loc_nodes
            local_dofs = (n -1) * num_dof_per_node .+ collect(1:num_dof_per_node)
            global_dofs = NODES[ELEMS[e].ChildNodes[n]].ChildDOFS
            GlobalExternalForceVector[global_dofs] += collect(ELEMS[e].ExternalForceVector[local_dofs])
        end
    end
    return GlobalExternalForceVector
end

function computeLocalElementStiffnessMatrices(ELEMS,NODES,D̃)
    num_elem = length(ELEMS)
    num_dof_per_node = length(NODES[1].ChildDOFS)
    for e = 1:num_elem
        nPts = 2
        side_id = 1
        eDegree = ELEMS[e].Degree
        num_loc_nodes = length(ELEMS[e].ChildNodes)
        ELEMS[e].StiffnessMatrix = NamedDimsArray{(:local_dof_id, :local_dof_id,)}(zeros(Float64, num_loc_nodes*num_dof_per_node, num_loc_nodes*num_dof_per_node))
        for qp_id = 1:size(ELEMS[e].Quadrature[1].QuadraturePoints,:local_qp_id)
            QP = ELEMS[e].Quadrature[side_id].QuadraturePoints[qp_id]
            # ξ,W = GaussQuadratureRule_2D(nPts) # QP.Coordinates
            for n1 = 1:num_loc_nodes
                n1_local_dofs = (n1-1) * num_dof_per_node .+ collect(1:num_dof_per_node)
                n1_global_dofs = NODES[ELEMS[e].ChildNodes[n1]].ChildDOFS   
                B̃₁ = StrainDisplacement_2D(∇LagrangeBasis_2D(eDegree[1],QP.Coordinates)[n1,:]) # ξ->StrainDisplacement_2D(∇LagrangeBasis_2D(eDegree[2],ξ)[n1,:]) #StrainDisplacement_2D(ELEMS[e].∇Nₐ[n1,:])
                for n2 = 1:num_loc_nodes
                    n2_local_dofs = (n2 -1) * num_dof_per_node .+ collect(1:num_dof_per_node)
                    n2_global_dofs = NODES[ELEMS[e].ChildNodes[n2]].ChildDOFS
                    B̃₂ = StrainDisplacement_2D(∇LagrangeBasis_2D(eDegree[2],QP.Coordinates)[n2,:]) # ξ->StrainDisplacement_2D(∇LagrangeBasis_2D(eDegree[1],ξ)[n2,:]) #StrainDisplacement_2D(ELEMS[e].∇Nₐ[n2,:])
                    ELEMS[e].StiffnessMatrix[n1_local_dofs,n2_local_dofs] += transpose(B̃₁) * D̃ * B̃₂ * QP.Weights
                end
            end
        end
    end
    return ELEMS
end

function assembleGlobalStiffnessMatrix(ELEMS,NODES)
    num_elem = length(ELEMS)
    num_nodes = length(NODES)
    num_dof_per_node = length(NODES[1].ChildDOFS)
    num_global_dofs = num_nodes * num_dof_per_node
    GlobalStiffnessMatrix = NamedDimsArray{(:global_dof_id, :global_dof_id,)}(zeros(Float64, num_global_dofs, num_global_dofs))
    for e = 1:num_elem
        num_loc_nodes = length(ELEMS[e].ChildNodes)
        for n1 = 1:num_loc_nodes
            n1_local_dofs = (n1 -1) * num_dof_per_node .+ collect(1:num_dof_per_node)
            n1_global_dofs = NODES[ELEMS[e].ChildNodes[n1]].ChildDOFS
            for n2 = 1:num_loc_nodes
                n2_local_dofs = (n2 -1) * num_dof_per_node .+ collect(1:num_dof_per_node)
                n2_global_dofs = NODES[ELEMS[e].ChildNodes[n2]].ChildDOFS
                GlobalStiffnessMatrix[n1_global_dofs, n2_global_dofs] += collect(ELEMS[e].StiffnessMatrix[n1_local_dofs,n2_local_dofs])
            end
        end
    end
    return GlobalStiffnessMatrix
end


function newton_raphson(ResidualFun, TangentFun, u₀, max_nl_step, max_newton_iter, ε) 
    nl_step = 0
    uₙ = u₀
    while nl_step < max_nl_step
        nl_step += 1
        newton_iter = 0
        R₀ = ResidualFun(uₙ)
        uᵢ = uₙ
        while newton_iter < max_newton_iter
            newton_iter += 1

            Rᵢ = ResidualFun(uᵢ)

            if LinearAlgebra.norm(Rᵢ) < ε*(LinearAlgebra.norm(R₀)) 
                break
            end

            K = TangentFun(uᵢ)

            Δu = collect(K) \ collect(Rᵢ)
            uᵢ += Δu  # Will need to replace with an "UpdateFun()"
        end
        uₙ = uᵢ
    end
    uₛ = uₙ
    return uₛ
end