module feQuadrature

import Roots
import DecFP
import ForwardDiff

include("feBasisFunctions.jl")
import .feBasisFunctions

export computeGaussQuadrature
export evaluateGaussQuadrature
function computeGaussQuadrature(nPt::Int)
    if nPt == 0
        ξᵢ = 0
        Wᵢ = 0
    else
        ξᵢ = Roots.find_zeros(ξ->feBasisFunctions.legendreBasis(nPt,DecFP.Dec128(ξ))[end],-1,1)
        Wᵢ = zeros(length(ξᵢ))
        f(ξ) = feBasisFunctions.legendreBasis(nPt,ξ)[end]
        for i = 1:length(Wᵢ)
            Wᵢ[i] = Float64(2 / ((1-DecFP.Dec128(ξᵢ[i])^2.) * DecFP.Dec128(ForwardDiff.derivative(f,DecFP.Dec128(ξᵢ[i]))^2.)))
        end
    end
    return ξᵢ , Wᵢ
end

function evaluateGaussQuadrature(fun,ξᵢ::Array{Float64,1},Wᵢ::Array{Float64,1})
    integral = fun(ξᵢ) .* transpose(Wᵢ)
    integral = sum(integral, dims=1)
    return integral
end
end