module feBasisFunctions
import LinearAlgebra
import ForwardDiff
import DecFP

function ⊗(f::Array{Float64},g::Array{Float64})
    return LinearAlgebra.kron(f,g)
end

function ⊗(f::Array{Any},g::Array{Any})
    return LinearAlgebra.kron(f,g)
end

function lagrangeBasis(degree::Any,variate::Any)
    node = range(-1,1,length=degree+1)
    L = Array{Any,1}(undef,degree+1)
    for j = 1:degree+1
        L[j] = 1.0
        for m = 1:degree+1
            if j != m
                L[j] *= ((variate-node[m])/(node[j]-node[m]))
            end
        end
    end
    return L
end

function lagrangeBasis(degree::Int64, variate::Float64)
    node = range(-1,1,length=degree+1)
    L = Array{Float64,1}(undef,degree+1)
    for j = 1:degree+1
        L[j] = 1.0
        for m = 1:degree+1
            if j != m
                L[j] *= ((variate-node[m])/(node[j]-node[m]))
            end
        end
    end
    return L
end

function legendreBasis(degree::Int64, ξ::Float64)
    L = zeros(degree+1)
    for n = 0:degree
        for k = 0:n
            L[n+1] = L[n+1] + binomial(n,k) * binomial(n+k,n) * ((ξ-1)/2)^k
        end
    end
    return L
end

function legendreBasis(degree::Int64, ξ::DecFP.Dec128)
    L = zeros(DecFP.Dec128, degree+1)
    for n = 0:degree
        for k = 0:n
            L[n+1] = L[n+1] + binomial(n,k) * binomial(n+k,n) * ((ξ-1)/2)^k
        end
    end
    return L
end

function legendreBasis(degree::Int64, ξ::Any)
    L = Array{Any,1}(undef, degree+1)
    for n = 0:degree
        L[n+1] = 0.
        for k = 0:n
            L[n+1] = L[n+1] + binomial(n,k) * binomial(n+k,n) * ((ξ-1)/2)^k
        end
    end
    return L
end

function bernsteinBasis(degree::Int64,ξ::Float64)
    scaleFactor = 2.0^(-degree) # Bi-unit domain scale factor from unit domain
    B = Array{Float64,1}(undef,degree+1)
    for i = 1:degree+1
        binCoeff = factorial(degree)/(factorial(i-1)*factorial(degree+1-i));
        B[i] = scaleFactor * binCoeff * (1-ξ)^(degree-(i-1))*(1+ξ)^(i-1);
    end
    return B
end

function bernsteinBasis(degree::Int64,ξ::Array{Float64,1})
    scaleFactor = 2.0^(-degree) # Bi-unit domain scale factor from unit domain
    B = Array{Float64,2}(undef,degree+1,length(ξ))
    for i = 1:degree+1
        binCoeff = factorial(degree)/(factorial(i-1)*factorial(degree+1-i));
        B[i,:] = scaleFactor * binCoeff * (1 .- ξ).^(degree-(i-1)) .* (1 .+ ξ).^(i-1);
    end
    return B
end

end