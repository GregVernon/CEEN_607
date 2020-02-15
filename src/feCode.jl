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
        ξ = [0.]
        W = [2.]
    elseif nPts == 2
        ξ = [-1/sqrt(3), 1/sqrt(3)]
        W = [1., 1.,]
    end
    return ξ,W
end

function GaussQuadratureRule_2D(nPts)
    W = Array{Array{Float64},1}(undef,2)
    ξ,W[1] = GaussQuadratureRule_1D(nPts)
    η,W[2] = GaussQuadratureRule_1D(nPts)

    nQP = nPts^2
    QP = NamedDimsArray{(:local_qp_id, :ℜ_2)}(zeros(Float64,nQP,2))
    QW = NamedDimsArray{(:local_qp_id, :qp_weight)}(zeros(Float64,nQP,1))
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
    QW = NamedDimsArray{(:local_qp_id, :qp_weight)}(zeros(Float64,nQP,1))
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