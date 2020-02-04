module feBasisFunctions
import LinearAlgebra
import ForwardDiff
import DecFP

include("feEnumerations.jl")
# include("feMesh.jl")
import .feEnumerations
# import ..feMesh

export ⊗

function ⊗(f::Array{Float64},g::Array{Float64})
    return LinearAlgebra.kron(f,g)
end

function ⊗(f::Array{Any},g::Array{Any})
    return LinearAlgebra.kron(f,g)
end

function constructBasis(::Val{feEnumerations.lagrange},elem_degree)
    num_dim = length(elem_degree)
    if     num_dim == 1
        u = ξ -> lagrangeBasis(elem_degree[1],ξ)
        return (ξ) -> u(ξ)
    elseif num_dim == 2
        u = ξ -> lagrangeBasis(elem_degree[1],ξ)
        v = η -> lagrangeBasis(elem_degree[2],η)
        return (ξ,η) -> u(ξ) ⊗ v(η)
    elseif num_dim == 3
        u = ξ -> lagrangeBasis(elem_degree[1],ξ)
        v = η -> lagrangeBasis(elem_degree[2],η)
        w = ζ -> lagrangeBasis(elem_degree[3],ζ)
        return (ξ,η,ζ) -> (u(ξ) ⊗ v(η)) ⊗ w(ζ)
    end
end

function constructBasis(::Val{feEnumerations.lagrange},elem_degree,loc_node_id)
    loc_node_indices = LocalNodeID_to_LocalNodeIndices(elem_degree,loc_node_id)
    println(loc_node_indices)
    num_dim = length(elem_degree)
    if     num_dim == 1
        u = ξ -> lagrangeBasis(elem_degree[1],loc_node_indices,ξ)
        return F(ξ) = u(ξ)
    elseif num_dim == 2
        u = ξ -> lagrangeBasis(elem_degree[1],loc_node_indices[1],ξ)
        v = η -> lagrangeBasis(elem_degree[2],loc_node_indices[2],η)
        return G(ξη) = u(ξη[1]) * v(ξη[2])
    elseif num_dim == 3
        u = ξ -> lagrangeBasis(elem_degree[1],loc_node_indices[1],ξ)
        v = η -> lagrangeBasis(elem_degree[2],loc_node_indices[2],η)
        w = ζ -> lagrangeBasis(elem_degree[3],loc_node_indices[3],ζ)
        return H(ξηζ) = (u(ξηζ[1]) * v(ξηζ[2])) * w(ξηζ[3])
    end
end

function constructBasis(::Val{feEnumerations.bernstein},elem_degree)
    num_dim = length(elem_degree)
    if     num_dim == 1
        u = ξ -> bernsteinBasis(elem_degree[1],ξ)
        return (ξ) -> u(ξ)
    elseif num_dim == 2
        u = ξ -> bernsteinBasis(elem_degree[1],ξ)
        v = η -> bernsteinBasis(elem_degree[2],η)
        return (ξ,η) -> u(ξ) ⊗ v(η)
    elseif num_dim == 3
        u = ξ -> bernsteinBasis(elem_degree[1],ξ)
        v = η -> bernsteinBasis(elem_degree[2],η)
        w = ζ -> bernsteinBasis(elem_degree[3],ζ)
        return (ξ,η,ζ) -> (u(ξ) ⊗ v(η)) ⊗ w(ζ)
    end
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

function lagrangeBasis(degree::Int64, index::Any, variate::Float64)
    node = range(-1,1,length=degree+1)
    L = 1.0
    for m = 1:degree+1
        if index != m
            L *= ((variate-node[m])/(node[index]-node[m]))
        end
    end
    return L
end

function lagrangeBasis(degree::Any, index::Any, variate::Any)
    node = range(-1,1,length=degree+1)
    L = 1.0
    for m = 1:degree+1
        if index != m
            L *= ((variate-node[m])/(node[index]-node[m]))
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

function LocalNodeID_to_LocalNodeIndices(elem_degree,loc_node_id)
    num_dim = length(elem_degree)
    if num_dim == 1
        return loc_node_id
    elseif num_dim == 2
        id = 0
        for ii = 1:elem_degree[1]+1
            for jj = 1:elem_degree[2]+1
                id+=1
                if id == loc_node_id
                    return [ii,jj]
                end
            end
        end
    elseif num_dim == 3
        id = 0
        for ii = 1:elem_degree[1]+1
            for jj = 1:elem_degree[2]+1
                for kk = 1:elem_degree[3]+1
                    id += 1
                    if id == loc_node_id
                        return [ii,jj,kk]
                    end
                end
            end
        end
    end
end

end