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

    ∫fun = 0.
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
    return ñ
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

            Δu = K \ Rᵢ
            uᵢ += Δu
        end
        uₙ = uᵢ
    end
    uₛ = uₙ
    return uₛ
end