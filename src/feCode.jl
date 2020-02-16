using NamedDims
import LinearAlgebra

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
    end
    return ξ,W
end

function GaussQuadratureRule_2D(nPts)
    W = Array{Array{Float64},1}(undef,2)
    ξ,W[1] = GaussQuadratureRule_1D(nPts)
    η,W[2] = GaussQuadratureRule_1D(nPts)

    nQP = nPts^2
    QP = NamedDimsArray{(:local_qp_id, :ℜ_2)}(zeros(Float64,nQP,2))
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
    QP = NamedDimsArray{(:local_qp_id, :ℜ_3)}(zeros(Float64,nQP,3))
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
    L = NamedDimsArray{(:local_node_id,)}(zeros(Float64,(deg+1)))
    if deg == 1
        L = [(1-ξ)/2, (1+ξ)/2]
        return L
    elseif deg == 2
        L = [(ξ^2 - ξ)/2, 1-ξ^2, (ξ^2 + ξ)/2]
        return L
    end
end

function LagrangeBasis_2D(deg,ξ)
    L = NamedDimsArray{(:local_node_id,)}(zeros(Float64,(deg+1)^2))
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
    L = NamedDimsArray{(:local_node_id,)}(zeros(Float64,(deg+1)^3))
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