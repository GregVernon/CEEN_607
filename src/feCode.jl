using NamedDims
import LinearAlgebra
import ForwardDiff

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
    QP = NamedDimsArray{(:local_qp_id, :ℙᴺ)}(zeros(Float64,nQP,2))
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

function LagrangeBasis_1D(deg, ξ)
    ξ = ξ[1]
    𝓝 = NamedDimsArray{(:local_node_id,)}(Array{Any,1}(undef,deg+1))
    if deg == 1
        𝓝 = [(1-ξ)/2, (1+ξ)/2]
        return 𝓝
    elseif deg == 2
        𝓝 = [(ξ^2 - ξ)/2, 1-ξ^2, (ξ^2 + ξ)/2]
        return 𝓝
    end
end

function LagrangeBasis_2D(deg,ξ)
    𝓝 = NamedDimsArray{(:local_node_id,)}(Array{Any,1}(undef,(deg+1)^2))
    𝓝₁ = LagrangeBasis_1D(deg,ξ[1])
    𝓝₂ = LagrangeBasis_1D(deg,ξ[2])
    n = 0
    for j = 1:deg+1
        for i = 1:deg+1
            n += 1
            𝓝[n] = 𝓝₁[i] * 𝓝₂[j]
        end
    end
    return 𝓝
end

function ∇LagrangeBasis_2D(deg,ξ)
    𝓝(ξ) = LagrangeBasis_2D(deg,ξ)
    ∇𝓝 = NamedDimsArray{(:local_node_id, :ℙᴺ,)}(ForwardDiff.jacobian(ξ->𝓝(ξ),ξ))
    return ∇𝓝
end

function map_ℙᴺ_to_ℝᴺ(𝓝, Xᵉ , ξ)
    x = NamedDimsArray{(:ℝᴺ,)}(zeros(Float64,size(ξ)))
    num_nodes = size(Xᵉ, :local_node_id)
    for n = 1:num_nodes
        x += 𝓝(ξ)[n] .* Xᵉ[n,:]
    end
    return x
end

function ∇map_ℙᴺ_to_ℝᴺ(∇𝓝, Xᵉ, ξ)
    num_nodes = size(Xᵉ,:local_node_id)
    num_par_dim = length(ξ)
    num_cart_dim = size(Xᵉ,:ℝᴺ)
    Jᵢⱼ = NamedDimsArray{(:ℝᴺ,:ℙᴺ)}(zeros(Float64,length(ξ),num_cart_dim))

    for n = 1:num_nodes
        for j = 1:num_cart_dim
            for i = 1:num_par_dim
                Jᵢⱼ[i,j] += ∇𝓝(ξ)[n,j] * Xᵉ[n,i] 
            end
        end
    end
    return Jᵢⱼ
end

function compute∇ₓ𝓝(∇𝓝, Jᵢⱼ)
    ∇ₓ𝓝 = ∇𝓝 * inv(Jᵢⱼ)
    return ∇ₓ𝓝
end

function StrainDisplacement_2D(∇𝓝ₐ)
    B = [∇𝓝ₐ[1]  0.0;
         0.0    ∇𝓝ₐ[2];
         ∇𝓝ₐ[2] ∇𝓝ₐ[1]]
end

function computeVirtualStrain_2D(∇𝓝, cₐ, ξ)
    num_nodes = length(cₐ)
    εᵥ = zeros(Float64,3)
    for a = 1:num_nodes
        εᵥ += StrainDisplacement_2D(∇𝓝[a,:]) * cₐ
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

    ñ = NamedDimsArray{(:ℙᴺ,)}(ñ[1:2])
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

function computeInternalForce(ũ, GEOM)
    GEOM.Elements = computeLocalInternalForceVector(ũ, GEOM.Elements, GEOM.Nodes)
    F_int = assembleGlobalInternalForceVector(GEOM.Elements, GEOM.Nodes)
end

function computeResidual(ũ, GEOM)
    F_External = computeExternalForce(GEOM)
    F_Internal = computeInternalForce(ũ, GEOM)
    Residual = F_External - F_Internal
    return Residual
end

function applyBoundaryConditions(R̃, K̃, ũ, GEOM)
    remove_dofs = Int[]
    num_dof_per_node = length(GEOM.Nodes[1].ChildDOFS)
    for ns_id = 1:length(GEOM.NodeSets)
        if isempty(GEOM.NodeSets[ns_id].BC_Type) == false
            for bc_id = 1:length(GEOM.NodeSets[ns_id].BC_Type)
                bc_type = GEOM.NodeSets[ns_id].BC_Type[bc_id]
                if Int(bc_type) == Int(feEnumerate.dirichlet)
                    for bc_dof = 1:length(GEOM.NodeSets[ns_id].BC_DOF)
                        constrained_local_dof = GEOM.NodeSets[ns_id].BC_DOF[bc_dof]
                        for n = 1:length(GEOM.NodeSets[ns_id].ChildNodes)
                            local_dofs = (n-1) * num_dof_per_node .+ collect(1:num_dof_per_node)
                            constrained_global_dofs = GEOM.Nodes[GEOM.NodeSets[ns_id].ChildNodes[n]].ChildDOFS[constrained_local_dof]
                            ũ[constrained_global_dofs] = GEOM.NodeSets[ns_id].BC_Value[bc_id]
                            append!(remove_dofs, constrained_global_dofs)
                        end
                    end
                end
            end
        end
    end
    keep_dofs = setdiff(1:length(R̃), remove_dofs)
    R̃ = R̃[keep_dofs]
    K̃ = K̃[keep_dofs, keep_dofs]
    return R̃, K̃, keep_dofs, ũ
end

function applyBoundaryConditions(ũ, GEOM)
    num_dof_per_node = length(GEOM.Nodes[1].ChildDOFS)
    for ns_id = 1:length(GEOM.NodeSets)
        if isempty(GEOM.NodeSets[ns_id].BC_Type) == false
            for bc_id = 1:length(GEOM.NodeSets[ns_id].BC_Type)
                bc_type = GEOM.NodeSets[ns_id].BC_Type[bc_id]
                α = GEOM.NodeSets[ns_id].BC_NL_Solve_Scale
                if Int(bc_type) == Int(feEnumerate.dirichlet)
                    for bc_dof = 1:length(GEOM.NodeSets[ns_id].BC_DOF)
                        constrained_local_dof = GEOM.NodeSets[ns_id].BC_DOF[bc_dof]
                        for n = 1:length(GEOM.NodeSets[ns_id].ChildNodes)
                            local_dofs = (n-1) * num_dof_per_node .+ collect(1:num_dof_per_node)
                            constrained_global_dofs = GEOM.Nodes[GEOM.NodeSets[ns_id].ChildNodes[n]].ChildDOFS[constrained_local_dof]
                            ũ[constrained_global_dofs] = α * GEOM.NodeSets[ns_id].BC_Value[bc_id]
                        end
                    end
                end
            end
        end
    end
    
    return ũ
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

    for e = 1:num_elem  # For each element...
        nPts = 2
        eDegree = ELEMS[e].Degree
        num_loc_nodes = length(ELEMS[e].ChildNodes)
        side_id = 1  # Integration occurs over the element interior
        for qp_id = 1:size(ELEMS[e].Quadrature[1].QuadraturePoints,:local_qp_id)  # For each quadrature point in the element... (Begin Gauss Quadrature)
            QP = ELEMS[e].Quadrature[side_id].QuadraturePoints[qp_id]  # (Let QP represent the quadrature point)
            for n = 1:num_loc_nodes  # For each node in the element...
                local_dofs = (n-1) * num_dof_per_node .+ collect(1:num_dof_per_node)  # Grab the element's local dof ids associated with the current node
                global_dofs = NODES[ELEMS[e].ChildNodes[n]].ChildDOFS  # Grab the global dof ids associated with the current node -- Should be a 1:1 mapping from local <-> global ids
                ∇𝓝 = ∇LagrangeBasis_2D(eDegree[1],QP.ℙ)
                ∇𝓝ₐ = ∇𝓝[n,:]
                B̃ = StrainDisplacement_2D(∇𝓝ₐ)
                ε = computeVirtualStrain_2D(∇𝓝, Δu[global_dofs], QP.ℙ)
                σ =  D̃ * ε
                ELEMS[e].InternalForceVector[local_dofs] += (transpose(B̃) * σ) * LinearAlgebra.det(QP.Jᵢⱼ) * QP.α * QP.𝒲
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
                                ELEMS[global_elem_id].ExternalForceVector[loc_dof_id] +=  f̃ * QP.Nₐ[loc_node_id] * QP.α * QP.𝒲
                            end
                        end
                    end
                end
            end
        end
    end

    # Evaluate forces on surface sets
    for ss_id = 1:length(SURFACESETS)  # Step through each surface set
        if isempty(SURFACESETS[ss_id].LC_Type) == false  # Check to see if at least one load is applied to this surface set
            for lc_id = 1:length(SURFACESETS[ss_id].LC_Type)  # For each load applied to this surface set...
                load_type = SURFACESETS[ss_id].LC_Type[lc_id]  # Get the type of the load
                if Int(load_type) == Int(feEnumerate.pressure)  # If the load is a pressure...
                    for loc_elem_id = 1:length(SURFACESETS[ss_id].ChildElements)  # For each element that has a surface (Element may appear multiple times)
                        global_elem_id = SURFACESETS[ss_id].ChildElements[loc_elem_id]  # Get the global element id
                        side_id = SURFACESETS[ss_id].ChildElements_LocalFace[loc_elem_id]  # Get the local side of the element that's included in the surface set
                        for qp_id = 1:length(ELEMS[global_elem_id].Quadrature[side_id].QuadraturePoints)  # For each quadrature point in the current side... (Begin Gauss Quadrature)
                            QP = ELEMS[global_elem_id].Quadrature[side_id].QuadraturePoints[qp_id]  # (Let QP represent the quadrature point)
                            f̃ = SURFACESETS[ss_id].LC_Magnitude[lc_id] * QP.ñ / LinearAlgebra.norm(QP.ñ)  # Compute the force vector associated with the pressure using the surface normal 
                            for loc_node_id = 1:length(ELEMS[global_elem_id].ChildNodes)  # For each node in the element...
                                loc_dof_id = (loc_node_id-1) * num_dof_per_node .+ collect(1:num_dof_per_node) # Grab the element's local dof ids associated with the current node
                                ELEMS[global_elem_id].ExternalForceVector[loc_dof_id] +=  f̃ * QP.𝓝[loc_node_id] * QP.α * QP.𝒲  # Add the GQ function evaluation to the local force vector  ### FIX ME -- Need change of coordinates
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
            # ξ,W = GaussQuadratureRule_2D(nPts) # QP.ℙ
            for n1 = 1:num_loc_nodes
                n1_local_dofs = (n1-1) * num_dof_per_node .+ collect(1:num_dof_per_node)
                n1_global_dofs = NODES[ELEMS[e].ChildNodes[n1]].ChildDOFS   
                B̃₁ = StrainDisplacement_2D(∇LagrangeBasis_2D(eDegree[1],QP.ℙ)[n1,:]) # ξ->StrainDisplacement_2D(∇LagrangeBasis_2D(eDegree[2],ξ)[n1,:]) #StrainDisplacement_2D(ELEMS[e].∇Nₐ[n1,:])
                for n2 = 1:num_loc_nodes
                    n2_local_dofs = (n2 -1) * num_dof_per_node .+ collect(1:num_dof_per_node)
                    n2_global_dofs = NODES[ELEMS[e].ChildNodes[n2]].ChildDOFS
                    B̃₂ = StrainDisplacement_2D(∇LagrangeBasis_2D(eDegree[2],QP.ℙ)[n2,:]) # ξ->StrainDisplacement_2D(∇LagrangeBasis_2D(eDegree[1],ξ)[n2,:]) #StrainDisplacement_2D(ELEMS[e].∇Nₐ[n2,:])
                    ELEMS[e].StiffnessMatrix[n1_local_dofs,n2_local_dofs] += transpose(B̃₁) * D̃ * B̃₂ * QP.𝒲
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

function InitialConditions(α, ũ, GEOM)
    for ns_id = 1:length(GEOM.NodeSets)
        if isempty(GEOM.NodeSets[ns_id].BC_Type) == false
            GEOM.NodeSets[ns_id].BC_NL_Solve_Scale = α
        end
    end

    u = applyBoundaryConditions(ũ, GEOM)
    R = computeResidual(u, GEOM)
    return R, u
end

function newton_raphson(GEOM, max_nl_step, max_newton_iter, ε) 
    # Initialize
    num_global_dofs = max(buildNodeGlobalDOFConnectivityArray(GEOM.Nodes)...)
    u₀ = zeros(Float64, num_global_dofs)
    nl_step = 0
    uₙ = u₀
    while nl_step < max_nl_step
        nl_step += 1
        newton_iter = 0
        R₀ = Inf64 #ResidualFun(uₙ)
        α = nl_step / max_nl_step;
        Rₙ, uₙ = InitialConditions(α, uₙ, GEOM)
        uᵢ = uₙ
        while newton_iter < max_newton_iter
            newton_iter += 1
            println("non-linear step: ", nl_step, " newton iteration: ", newton_iter)
            Rᵢ = computeResidual(uᵢ,GEOM)
            println(Rᵢ)
            if LinearAlgebra.norm(Rᵢ) < ε 
                break
            end
            println("  residual: ", LinearAlgebra.norm(Rᵢ))

            K = assembleGlobalStiffnessMatrix(GEOM.Elements, GEOM.Nodes);
            Rᵢ, K, keep_dofs, uᵢ = applyBoundaryConditions(Rᵢ, K, uᵢ, GEOM)

            Δu = collect(K) \ collect(Rᵢ)
            uᵢ[keep_dofs] += Δu  # Will need to replace with an "UpdateFun()"
        end
        uₙ = uᵢ
    end
    uₛ = uₙ
    return uₛ
end