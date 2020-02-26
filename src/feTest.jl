using Test
import Base.≈
≈(a::Float64,b::Float64) = isapprox(a::Float64,b::Float64;atol=10*eps(),rtol=10*eps())

include("feCode.jl")
include("feElement.jl")
include("feUtility.jl")

# Test 1D Guass Quadrature
@testset "1D Guass Quadrature" begin
    @testset "1 Point Rule" begin
        nPts = 1
        @testset "Degree = 0" begin
            @test GaussQuadrature_1D(x->(x^0), nPts) ≈ 2.0
        end
        @testset "Degree = 1" begin
            @test GaussQuadrature_1D(x->(x^1), nPts) ≈ 0.0
        end
    end

    @testset "1 Point Rule" begin
        nPts = 2
        @testset "Degree = 0" begin
            @test GaussQuadrature_1D(x->(x^0), nPts) ≈ 2.0
        end
        @testset "Degree = 1" begin
            @test GaussQuadrature_1D(x->(x^1), nPts) ≈ 0.0
        end
        @testset "Degree = 2" begin
            @test GaussQuadrature_1D(x->(x^2), nPts) ≈ 2.0/3.0
        end
        @testset "Degree = 3" begin
            @test GaussQuadrature_1D(x->(x^3), nPts) ≈ 0.
        end
    end
end

# Test 2D Guass Quadrature
@testset "2D Guass Quadrature" begin
    @testset "1x1 Point Rule" begin
        nPts = 1
        @testset "Degree = 0" begin
            @test GaussQuadrature_2D((xy)->(xy[1]^0 * xy[2]^0), nPts) ≈ 4.0
        end
        @testset "Degree = 1" begin
            exInt = [4., 0., 0., 0.]
            n = 0
            for j = 0:1
                for i = 0:1
                    n+=1
                    @test GaussQuadrature_2D((xy)->(xy[1]^i * xy[2]^j), nPts) ≈ exInt[n]
                end
            end 
        end
    end

    @testset "2x2 Point Rule" begin
        nPts = 2
        @testset "Degree = 0" begin
            @test GaussQuadrature_2D((xy)->(xy[1]^0 * xy[2]^0), nPts) ≈ 4.0
        end
        @testset "Degree = 1" begin
            exInt = [4., 0., 0., 0.]
            n = 0
            for j = 0:1
                for i = 0:1
                    n+=1
                    @test GaussQuadrature_2D((xy)->(xy[1]^i * xy[2]^j), nPts) ≈ exInt[n]
                end
            end 
        end
        @testset "Degree = 2" begin
            exInt = [4., 0., 4/3, 0., 0., 0., 4/3, 0., 4/9]
            n = 0
            for j = 0:2
                for i = 0:2
                    n+=1
                    @test GaussQuadrature_2D((xy)->(xy[1]^i * xy[2]^j), nPts) ≈ exInt[n]
                end
            end
        end
        @testset "Degree = 3" begin
            exInt = [4., 0., 4/3, 0., 0., 0., 0., 0., 4/3, 0., 4/9, 0., 0., 0., 0., 0.]
            n = 0
            for j = 0:3
                for i = 0:3
                    n+=1
                    @test GaussQuadrature_2D((xy)->(xy[1]^i * xy[2]^j), nPts) ≈ exInt[n]
                end
            end
        end
    end
end

# Test 3D Guass Quadrature
@testset "3D Guass Quadrature" begin
    @testset "1x1x1 Point Rule" begin
        nPts = 1
        @testset "Degree = 0" begin
            @test GaussQuadrature_3D((xyz)->(xyz[1]^0 * xyz[2]^0 * xyz[3]^0), nPts) ≈ 8.0
        end
        @testset "Degree = 1" begin
            exInt = [8.0, 0, 0, 0, 0, 0, 0, 0, 0]
            n = 0
            for k = 0:1
                for j = 0:1
                    for i = 0:1
                        n += 1
                        @test GaussQuadrature_3D((xyz)->(xyz[1]^i * xyz[2]^j * xyz[3]^k), nPts) ≈ exInt[n] # 0
                    end
                end
            end
        end
    end

    @testset "2x2x2 Point Rule" begin
        nPts = 2
        @testset "Degree = 0" begin
            @test GaussQuadrature_3D((xyz)->(xyz[1]^0 * xyz[2]^0 * xyz[3]^0), nPts) ≈ 8.0
        end
        @testset "Degree = 1" begin
            exInt = [8.0, 0, 0, 0, 0, 0, 0, 0, 0]
            n = 0
            for k = 0:1
                for j = 0:1
                    for i = 0:1
                        n += 1
                        @test GaussQuadrature_3D((xyz)->(xyz[1]^i * xyz[2]^j * xyz[3]^k), nPts) ≈ exInt[n] # 0
                    end
                end
            end
        end
        @testset "Degree = 2" begin
            exInt = [8.0, 0, 8/3, 0, 0, 0, 8/3, 0, 8/9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8/3, 0, 8/9, 0, 0, 0, 8/9, 0, 8/27]
            n = 0
            for k = 0:2
                for j = 0:2
                    for i = 0:2
                        n+=1
                        @test GaussQuadrature_3D((xyz)->(xyz[1]^i * xyz[2]^j * xyz[3]^k), nPts) ≈ exInt[n]
                    end
                end
            end            
        end
        @testset "Degree = 3" begin
            exInt = zeros(Float64,81)
            exInt[1] = 8.
            exInt[3] = 8/3
            exInt[9] = 8/3
            exInt[11] = 8/9
            exInt[33] = 8/3
            exInt[35] = 8/9
            exInt[41] = 8/9
            exInt[43] = 8/27
            n = 0
            for k = 0:3
                for j = 0:3
                    for i = 0:3
                        n+=1
                        @test GaussQuadrature_3D((xyz)->(xyz[1]^i * xyz[2]^j * xyz[3]^k), nPts) ≈ exInt[n]
                    end
                end
            end   
        end
    end
end

# Test 1D Lagrange Basis Function
@testset "1D Lagrange Basis Function" begin
    @testset "Degree = 1" begin
        degree = 1
        numNodes = degree + 1
        nodeCoords = [-1., 1.]
        for n = 1:numNodes
            @test LagrangeBasis_1D(degree, nodeCoords[n])[n] ≈ 1.0
        end
    end

    @testset "Degree = 2" begin
        degree = 2
        numNodes = degree + 1
        nodeCoords = [-1., 0., 1.]
        for n = 1:numNodes
            @test LagrangeBasis_1D(degree, nodeCoords[n])[n] ≈ 1.0
        end
    end
end

# Test 2D Lagrange Basis Function
@testset "2D Lagrange Basis Function" begin
    @testset "Degree = 1" begin
        degree = 1
        numNodes = (degree + 1) ^ 2
        nodeCoords = buildLocalNodeCoordinates_2D(degree)
        for n = 1:numNodes
            @test LagrangeBasis_2D(degree, nodeCoords[n,:])[n] ≈ 1.0
        end
    end

    @testset "Degree = 2" begin
        degree = 2
        numNodes = (degree + 1)^2
        nodeCoords = buildLocalNodeCoordinates_2D(degree)
        for n = 1:numNodes
            @test LagrangeBasis_2D(degree, nodeCoords[n,:])[n] ≈ 1.0
        end
    end
end

# Test 3D Lagrange Basis Function
@testset "3D Lagrange Basis Function" begin
    @testset "Degree = 1" begin
        degree = 1
        numNodes = (degree + 1) ^ 3
        nodeCoords = buildLocalNodeCoordinates_3D(degree)
        for n = 1:numNodes
            @test LagrangeBasis_3D(degree, nodeCoords[n,:])[n] ≈ 1.0
        end
    end

    @testset "Degree = 2" begin
        degree = 2
        numNodes = (degree + 1)^3
        nodeCoords = buildLocalNodeCoordinates_3D(degree)
        for n = 1:numNodes
            @test LagrangeBasis_3D(degree, nodeCoords[n,:])[n] ≈ 1.0
        end
    end
end

@testset "1D Lagrange Basis Derivatives" begin
    @testset "Degree = 1" begin
        degree = 1
        @test ∇LagrangeBasis_1D(degree,[-1.0]) ≈ [-0.5, 0.5]
        @test ∇LagrangeBasis_1D(degree,[+1.0]) ≈ [-0.5, 0.5]
    end

    @testset "Degree = 2" begin
        degree = 2
        @test ∇LagrangeBasis_1D(degree,[-1.0]) ≈ [-3/2, 2, -1/2]
        @test ∇LagrangeBasis_1D(degree,[0.0])  ≈ [-1/2, 0, 1/2]
        @test ∇LagrangeBasis_1D(degree,[+1.0]) ≈ [1/2, -2, 3/2]
    end
end

@testset "2D Lagrange Basis Derivatives" begin
    @testset "Degree = 1" begin
        degree = 1
        @test ∇LagrangeBasis_2D(degree,[-1.0, -1.0]) ≈ [-0.5 -0.5; 
                                                         0.5  0.0; 
                                                         0.0  0.5; 
                                                         0.0  0.0]
        @test ∇LagrangeBasis_2D(degree,[+1.0, +1.0]) ≈ [ 0.0  0.0; 
                                                         0.0 -0.5; 
                                                        -0.5  0.0; 
                                                         0.5  0.5]
    end

    @testset "Degree = 2" begin
        degree = 2
        @test ∇LagrangeBasis_2D(degree,[-1.0, -1.0]) ≈ [-3/2 -3/2;
                                                         2.    0.;
                                                        -1/2   0.;
                                                         0.    2.;
                                                         0.    0.;
                                                         0.    0.;
                                                         0.  -1/2;
                                                         0.    0.;
                                                         0.    0.]
        @test ∇LagrangeBasis_2D(degree,[0.0, 0.0])   ≈ [ 0.    0.;
                                                         0.  -1/2;
                                                         0     0.;
                                                        -1/2.  0.;
                                                         0.    0.;
                                                         1/2.  0.;
                                                         0.    0.;
                                                         0.   1/2;
                                                         0.    0.]
        @test ∇LagrangeBasis_2D(degree,[+1.0, +1.0]) ≈ [ 0.    0.;
                                                         0.    0.;
                                                         0.   1/2;
                                                         0.    0.;
                                                         0.    0.;
                                                         0.   -2.;
                                                         1/2   0.;
                                                        -2.    0.;
                                                         3/2  3/2]
    end
end

@testset "3D Lagrange Basis Derivatives" begin
    @testset "Degree = 1" begin
        degree = 1
        @test ∇LagrangeBasis_3D(degree,[-1.0, -1.0, -1.0]) ≈ [-1/2 -1/2 -1/2;
                                                               1/2   0.   0.;
                                                               0.   1/2   0.;
                                                               0.    0.   0.;
                                                               0.    0.  1/2; 
                                                               0.    0.   0.;
                                                               0.    0.   0.;
                                                               0.    0.   0.]

        @test ∇LagrangeBasis_3D(degree,[+1.0, +1.0, +1.0]) ≈ [ 0.   0.   0.; 
                                                               0.   0.   0.;
                                                               0.   0.   0.;
                                                               0.   0. -1/2;
                                                               0.   0.   0.;
                                                               0. -1/2   0.;
                                                             -1/2   0.   0.;
                                                              1/2  1/2  1/2]
                                                         
    end

    @testset "Degree = 2" begin
        degree = 2
        exJac = zeros(Float64,27,3)
        exJac[1,:] = [-3/2 -3/2 -3/2]
        exJac[2,1] = 2.
        exJac[3,1] = -1/2
        exJac[4,2] = 2.
        exJac[7,2] = -1/2
        exJac[10,3] = 2.
        exJac[19,3] = -1/2
        @test ∇LagrangeBasis_3D(degree,[-1.0, -1.0, -1.0]) ≈ exJac
        exJac = zeros(Float64,27,3)
        exJac[5,3]  = -1/2
        exJac[11,2] = -1/2
        exJac[13,1] = -1/2
        exJac[15,1] = +1/2
        exJac[17,2] = +1/2
        exJac[23,3] = +1/2
        @test ∇LagrangeBasis_3D(degree,[0.0, 0.0, 0.0])   ≈ exJac
        exJac = zeros(Float64,27,3)
        exJac[9,3]  = 1/2 
        exJac[18,3] = -2
        exJac[21,2] = 1/2
        exJac[24,2] = -2
        exJac[25,1] = 1/2
        exJac[26,1] = -2
        exJac[27,:] = [3/2 3/2 3/2]
        @test ∇LagrangeBasis_3D(degree,[+1.0, +1.0, +1.0]) ≈ exJac
    end
end

@testset "Geometric Mapping" begin
    @testset "Dimension = 1" begin
        @testset "Degree = 1" begin
            deg = 1
            Nₐ = ξ->LagrangeBasis_1D(deg,ξ)
            xₐ = buildLocalNodeCoordinates_1D(deg)
            ξ = rand(1)
            @test all(computeGeometricMapping(Nₐ, xₐ, ξ) .≈ ξ)
            @test all(computeGeometricMapping(Nₐ, 2*xₐ, ξ) .≈ 2*ξ)
            @test all(computeGeometricMapping(Nₐ, xₐ.+2, ξ) .≈ ξ .+ 2)
        end

        @testset "Degree = 2" begin
            deg = 2
            Nₐ = ξ->LagrangeBasis_1D(deg,ξ)
            xₐ = buildLocalNodeCoordinates_1D(deg)
            ξ = rand(1)
            @test all(computeGeometricMapping(Nₐ, xₐ, ξ) .≈ ξ)
            @test all(computeGeometricMapping(Nₐ, 2*xₐ, ξ) .≈ 2*ξ)
            @test all(computeGeometricMapping(Nₐ, xₐ.+2, ξ) .≈ ξ .+ 2)
        end
    end

    @testset "Dimension = 2" begin
        @testset "Degree = 1" begin
            deg = 1
            Nₐ = ξ->LagrangeBasis_2D(deg,ξ)
            xₐ = buildLocalNodeCoordinates_2D(deg)
            ξ = rand(2)
            @test all(computeGeometricMapping(Nₐ, xₐ, ξ) .≈ ξ)
            @test all(computeGeometricMapping(Nₐ, 2*xₐ, ξ) .≈ 2*ξ)
            @test all(computeGeometricMapping(Nₐ, xₐ.+2, ξ) .≈ ξ .+ 2)
        end

        @testset "Degree = 2" begin
            deg = 2
            Nₐ = ξ->LagrangeBasis_2D(deg,ξ)
            xₐ = buildLocalNodeCoordinates_2D(deg)
            ξ = rand(2)
            @test all(computeGeometricMapping(Nₐ, xₐ, ξ) .≈ ξ)
            @test all(computeGeometricMapping(Nₐ, 2*xₐ, ξ) .≈ 2*ξ)
            @test all(computeGeometricMapping(Nₐ, xₐ.+2, ξ) .≈ ξ .+ 2)
        end
    end

    @testset "Dimension = 3" begin
        @testset "Degree = 1" begin
            deg = 1
            Nₐ = ξ->LagrangeBasis_3D(deg,ξ)
            xₐ = buildLocalNodeCoordinates_3D(deg)
            ξ = rand(3)
            @test all(computeGeometricMapping(Nₐ, xₐ, ξ) .≈ ξ)
            @test all(computeGeometricMapping(Nₐ, 2*xₐ, ξ) .≈ 2*ξ)
            @test all(computeGeometricMapping(Nₐ, xₐ.+2, ξ) .≈ ξ .+ 2)
        end

        @testset "Degree = 2" begin
            deg = 2
            Nₐ = ξ->LagrangeBasis_3D(deg,ξ)
            xₐ = buildLocalNodeCoordinates_3D(deg)
            ξ = rand(3)
            @test all(computeGeometricMapping(Nₐ, xₐ, ξ) .≈ ξ)
            @test all(computeGeometricMapping(Nₐ, 2*xₐ, ξ) .≈ 2*ξ)
            @test all(computeGeometricMapping(Nₐ, xₐ.+2, ξ) .≈ ξ .+ 2)
        end
    end
end

@testset "Jacobian of Geometric Mapping" begin
    @testset "Dimension = 1" begin
        @testset "Degree = 1" begin
            degree = 1
            ∇Nₐ = ξ->∇LagrangeBasis_1D(degree,ξ)
            xₐ = buildLocalNodeCoordinates_1D(degree)
            for n = 1:10
                ξ = 2*rand(1) .- 1
                @test compute∇GeometricMapping(∇Nₐ, xₐ, ξ) ≈ NamedDimsArray{(:ℝᴺ,:ℙᴺ,)}(LinearAlgebra.I(1))
            end
        end

        @testset "Degree = 2" begin
            degree = 2
            ∇Nₐ = ξ->∇LagrangeBasis_1D(degree,ξ)

            xₐ = buildLocalNodeCoordinates_1D(degree)
            for n = 1:10
                ξ = 2*rand(1) .- 1
                @test compute∇GeometricMapping(∇Nₐ, xₐ, ξ) ≈ NamedDimsArray{(:ℝᴺ,:ℙᴺ,)}(LinearAlgebra.I(1))
            end

            xₐ = 2 * buildLocalNodeCoordinates_1D(degree)
            for n = 1:10
                ξ = 2*(2*rand(1) .- 1)
                @test compute∇GeometricMapping(∇Nₐ, xₐ, ξ) ≈ 2 * NamedDimsArray{(:ℝᴺ,:ℙᴺ,)}(LinearAlgebra.I(1))
            end

            xₐ = 2 .+ buildLocalNodeCoordinates_1D(degree)
            for n = 1:10
                ξ = 2 .+ (2*rand(1) .- 1)
                @test compute∇GeometricMapping(∇Nₐ, xₐ, ξ) ≈ NamedDimsArray{(:ℝᴺ,:ℙᴺ,)}(LinearAlgebra.I(1))
            end

            xₐ = 2 .+ (2*buildLocalNodeCoordinates_1D(degree))
            for n = 1:10
                ξ = 2 .+ 2*(2*rand(1) .- 1)
                @test compute∇GeometricMapping(∇Nₐ, xₐ, ξ) ≈ 2 * NamedDimsArray{(:ℝᴺ,:ℙᴺ,)}(LinearAlgebra.I(1))
            end
        end
    end

    @testset "Dimension = 2" begin
        @testset "Degree = 1" begin
            degree = 1
            ∇Nₐ = ξ->∇LagrangeBasis_2D(degree,ξ)
            
            xₐ = buildLocalNodeCoordinates_2D(degree)
            for n = 1:10
                ξ = 2*rand(2) .- 1
                @test compute∇GeometricMapping(∇Nₐ, xₐ, ξ) ≈ NamedDimsArray{(:ℝᴺ,:ℙᴺ,)}(LinearAlgebra.I(2))
            end
            
            xₐ = 2 .* buildLocalNodeCoordinates_2D(degree)
            for n = 1:10
                ξ = 2*(2*rand(2) .- 1)
                @test compute∇GeometricMapping(∇Nₐ, xₐ, ξ) ≈ 2 * NamedDimsArray{(:ℝᴺ,:ℙᴺ,)}(LinearAlgebra.I(2))
            end
            
            xₐ = 2 .+ buildLocalNodeCoordinates_2D(degree)
            for n = 1:10
                ξ = 2 .+ (2*rand(2) .- 1)
                @test compute∇GeometricMapping(∇Nₐ, xₐ, ξ) ≈ NamedDimsArray{(:ℝᴺ,:ℙᴺ,)}(LinearAlgebra.I(2))
            end
            
            xₐ = 2 .+ (2*buildLocalNodeCoordinates_2D(degree))
            for n = 1:10
                ξ = 2 .+ 2*(2*rand(2) .- 1)
                @test compute∇GeometricMapping(∇Nₐ, xₐ, ξ) ≈ 2 * NamedDimsArray{(:ℝᴺ,:ℙᴺ,)}(LinearAlgebra.I(2))
            end
        end

        @testset "Degree = 2" begin
            degree = 2
            ∇Nₐ = ξ->∇LagrangeBasis_2D(degree,ξ)
            
            xₐ = buildLocalNodeCoordinates_2D(degree)
            for n = 1:10
                ξ = 2*rand(2) .- 1
                @test compute∇GeometricMapping(∇Nₐ, xₐ, ξ) ≈ NamedDimsArray{(:ℝᴺ,:ℙᴺ,)}(LinearAlgebra.I(2))
            end
            
            xₐ = 2 .* buildLocalNodeCoordinates_2D(degree)
            for n = 1:10
                ξ = 2*(2*rand(2) .- 1)
                @test compute∇GeometricMapping(∇Nₐ, xₐ, ξ) ≈ 2 * NamedDimsArray{(:ℝᴺ,:ℙᴺ,)}(LinearAlgebra.I(2))
            end
            
            xₐ = 2 .+ buildLocalNodeCoordinates_2D(degree)
            for n = 1:10
                ξ = 2 .+ (2*rand(2) .- 1)
                @test compute∇GeometricMapping(∇Nₐ, xₐ, ξ) ≈ NamedDimsArray{(:ℝᴺ,:ℙᴺ,)}(LinearAlgebra.I(2))
            end
            
            xₐ = 2 .+ (2*buildLocalNodeCoordinates_2D(degree))
            for n = 1:10
                ξ = 2 .+ 2*(2*rand(2) .- 1)
                @test compute∇GeometricMapping(∇Nₐ, xₐ, ξ) ≈ 2 * NamedDimsArray{(:ℝᴺ,:ℙᴺ,)}(LinearAlgebra.I(2))
            end
        end
    end

    @testset "Dimension = 3" begin
        @testset "Degree = 1" begin
            degree = 1
            ∇Nₐ = ξ->∇LagrangeBasis_3D(degree,ξ)
            
            xₐ = buildLocalNodeCoordinates_3D(degree)
            for n = 1:10
                ξ = 2*rand(3) .- 1
                @test compute∇GeometricMapping(∇Nₐ, xₐ, ξ) ≈ NamedDimsArray{(:ℝᴺ,:ℙᴺ,)}(LinearAlgebra.I(3))
            end

            xₐ = 2. * buildLocalNodeCoordinates_3D(degree)
            for n = 1:10
                ξ = 2*(2*rand(3) .- 1)
                @test compute∇GeometricMapping(∇Nₐ, xₐ, ξ) ≈ 2 * NamedDimsArray{(:ℝᴺ,:ℙᴺ,)}(LinearAlgebra.I(3))
            end

            xₐ = 2. .+ buildLocalNodeCoordinates_3D(degree)
            for n = 1:10
                ξ = 2 .+ (2*rand(3) .- 1)
                @test compute∇GeometricMapping(∇Nₐ, xₐ, ξ) ≈ NamedDimsArray{(:ℝᴺ,:ℙᴺ,)}(LinearAlgebra.I(3))
            end

            xₐ = 2. .+ (2*buildLocalNodeCoordinates_3D(degree))
            for n = 1:10
                ξ = 2 .+ 2*((2*rand(3) .- 1))
                @test compute∇GeometricMapping(∇Nₐ, xₐ, ξ) ≈ 2 * NamedDimsArray{(:ℝᴺ,:ℙᴺ,)}(LinearAlgebra.I(3))
            end
        end

        @testset "Degree = 2" begin
            degree = 2
            ∇Nₐ = ξ->∇LagrangeBasis_3D(degree,ξ)

            xₐ = buildLocalNodeCoordinates_3D(degree)
            for n = 1:10
                ξ = 2*rand(3) .- 1
                @test compute∇GeometricMapping(∇Nₐ, xₐ, ξ) ≈ NamedDimsArray{(:ℝᴺ,:ℙᴺ,)}(LinearAlgebra.I(3))
            end

            xₐ = 2. * buildLocalNodeCoordinates_3D(degree)
            for n = 1:10
                ξ = 2*(2*rand(3) .- 1)
                @test compute∇GeometricMapping(∇Nₐ, xₐ, ξ) ≈ 2 * NamedDimsArray{(:ℝᴺ,:ℙᴺ,)}(LinearAlgebra.I(3))
            end

            xₐ = 2. .+ buildLocalNodeCoordinates_3D(degree)
            for n = 1:10
                ξ = 2 .+ (2*rand(3) .- 1)
                @test compute∇GeometricMapping(∇Nₐ, xₐ, ξ) ≈ NamedDimsArray{(:ℝᴺ,:ℙᴺ,)}(LinearAlgebra.I(3))
            end

            xₐ = 2. .+ (2*buildLocalNodeCoordinates_3D(degree))
            for n = 1:10
                ξ = 2 .+ 2*((2*rand(3) .- 1))
                @test compute∇GeometricMapping(∇Nₐ, xₐ, ξ) ≈ 2 * NamedDimsArray{(:ℝᴺ,:ℙᴺ,)}(LinearAlgebra.I(3))
            end
        end
    end
end

@testset "Boundary Normal" begin
    @testset "Dimension = 2" begin
        @testset "Degree = 1" begin
            degree = 1
            xₐ = buildLocalNodeCoordinates_2D(degree)
            ∇Nₐ = ξ->∇LagrangeBasis_2D(degree,ξ)
            Jᵢⱼ = ξ->compute∇GeometricMapping(∇Nₐ, xₐ, ξ)
            for sideID = 1:4
                if     sideID == 1
                    ξ = [-1.0, 0.0]
                    ñ = [-1.0, 0.0, 0.0]
                elseif sideID == 2
                    ξ = [+1.0, 0.0]
                    ñ = [+1.0, 0.0, 0.0]
                elseif sideID == 3
                    ξ = [0.0, -1.0]
                    ñ = [0.0, -1.0, 0.0]
                else   sideID == 4
                    ξ = [0.0, +1.0]
                    ñ = [0.0, +1.0, 0.0]
                end
                n = computeBoundaryNormals(Jᵢⱼ(ξ), sideID)
                @test n ≈ ñ
            end

            A = affineTransformationMatrix_2D(rotate=pi/4)
            xₐ = buildLocalNodeCoordinates_2D(degree)
            for n = 1:size(xₐ,:local_node_id)
                xₐ[n,:] = (A * [xₐ[n,:]...,0.0])[1:2]
            end
            ∇Nₐ = ξ->∇LagrangeBasis_2D(degree,ξ)
            Jᵢⱼ = ξ->compute∇GeometricMapping(∇Nₐ, xₐ, ξ)
            for sideID = 1:4
                if     sideID == 1
                    ξ = [-1.0, 0.0]
                    ñ = [-√2/2, -√2/2, 0.0]
                elseif sideID == 2
                    ξ = [+1.0, 0.0]
                    ñ = [+√2/2, +√2/2, 0.0]
                elseif sideID == 3
                    ξ = [0.0, -1.0]
                    ñ = [+√2/2, -√2/2, 0.0]
                else   sideID == 4
                    ξ = [0.0, +1.0]
                    ñ = [-√2/2, +√2/2, 0.0]
                end
                n = computeBoundaryNormals(Jᵢⱼ(ξ), sideID)
                @test n ≈ ñ
            end
        end

        @testset "Degree = 2" begin
            degree = 2
            xₐ = buildLocalNodeCoordinates_2D(degree)
            ∇Nₐ = ξ->∇LagrangeBasis_2D(degree,ξ)
            Jᵢⱼ = ξ->compute∇GeometricMapping(∇Nₐ, xₐ, ξ)
            for sideID = 1:4
                if     sideID == 1
                    ξ = [-1.0, 0.0]
                    ñ = [-1.0, 0.0, 0.0]
                elseif sideID == 2
                    ξ = [+1.0, 0.0]
                    ñ = [+1.0, 0.0, 0.0]
                elseif sideID == 3
                    ξ = [0.0, -1.0]
                    ñ = [0.0, -1.0, 0.0]
                else   sideID == 4
                    ξ = [0.0, +1.0]
                    ñ = [0.0, +1.0, 0.0]
                end
                n = computeBoundaryNormals(Jᵢⱼ(ξ), sideID)
                @test n ≈ ñ
            end

            A = affineTransformationMatrix_2D(rotate=pi/4)
            xₐ = buildLocalNodeCoordinates_2D(degree)
            for n = 1:size(xₐ,:local_node_id)
                xₐ[n,:] = (A * [xₐ[n,:]...,0.0])[1:2]
            end
            ∇Nₐ = ξ->∇LagrangeBasis_2D(degree,ξ)
            Jᵢⱼ = ξ->compute∇GeometricMapping(∇Nₐ, xₐ, ξ)
            for sideID = 1:4
                if     sideID == 1
                    ξ = [-1.0, 0.0]
                    ñ = [-√2/2, -√2/2, 0.0]
                elseif sideID == 2
                    ξ = [+1.0, 0.0]
                    ñ = [+√2/2, +√2/2, 0.0]
                elseif sideID == 3
                    ξ = [0.0, -1.0]
                    ñ = [+√2/2, -√2/2, 0.0]
                else   sideID == 4
                    ξ = [0.0, +1.0]
                    ñ = [-√2/2, +√2/2, 0.0]
                end
                n = computeBoundaryNormals(Jᵢⱼ(ξ), sideID)
                @test n ≈ ñ
            end
        end
    end
end

@testset "Integral Scaling" begin
    @testset "2D Element" begin
        @testset "Degree = 1" begin
            degree = 1
            xₐ = buildLocalNodeCoordinates_2D(degree)
            ∇Nₐ = ξ->∇LagrangeBasis_2D(degree,ξ)
            Jᵢⱼ = ξ->compute∇GeometricMapping(∇Nₐ, xₐ, ξ)
            @testset "Body" begin
                sideID = 0
                ξ = [0.0, 0.0]
                IntScale = computeIntegralScaling_2D(Jᵢⱼ(ξ), sideID)
                @test IntScale ≈ 1.0
            end
            @testset "Sides" begin
                for sideID = 1:4
                    if     sideID == 1
                        ξ = [-1.0, 0.0]
                        IntScale = computeIntegralScaling_2D(Jᵢⱼ(ξ), sideID)
                    elseif sideID == 2
                        ξ = [+1.0, 0.0]
                        IntScale = computeIntegralScaling_2D(Jᵢⱼ(ξ), sideID)
                    elseif sideID == 3
                        ξ = [0.0, -1.0]
                        IntScale = computeIntegralScaling_2D(Jᵢⱼ(ξ), sideID)
                    elseif sideID == 4
                        ξ = [0.0, +1.0]
                        IntScale = computeIntegralScaling_2D(Jᵢⱼ(ξ), sideID)
                    end
                    @test IntScale ≈ 1.0
                end
            end
        end
    end
end

@testset "Change of Coordinates" begin
    @testset "Dimension = 1" begin
        @testset "Degree = 1" begin
            degree = 1
            ∇Nₐ = ξ->∇LagrangeBasis_1D(degree,ξ)
            
            xₐ = buildLocalNodeCoordinates_1D(degree)
            Jᵢⱼ = ξ->compute∇GeometricMapping(∇Nₐ, xₐ, ξ)
            ξ = [0.0]
            ∇ₓNₐ = compute∇ₓNₐ(∇Nₐ(ξ), Jᵢⱼ(ξ))
            @test ∇ₓNₐ ≈ ∇Nₐ(ξ)
        end
        @testset "Degree = 2" begin
            degree = 2
            ∇Nₐ = ξ->∇LagrangeBasis_1D(degree,ξ)

            xₐ = buildLocalNodeCoordinates_1D(degree)
            Jᵢⱼ = ξ->compute∇GeometricMapping(∇Nₐ, xₐ, ξ)
            ξ = [0.0]
            ∇ₓNₐ = compute∇ₓNₐ(∇Nₐ(ξ), Jᵢⱼ(ξ))
            @test ∇ₓNₐ ≈ ∇Nₐ(ξ)
        end
    end

    @testset "Dimension = 2" begin
        @testset "Degree = 1" begin
            degree = 1
            ∇Nₐ = ξ->∇LagrangeBasis_2D(degree,ξ)
            
            # Identity
            xₐ = buildLocalNodeCoordinates_2D(degree)
            Jᵢⱼ = ξ->compute∇GeometricMapping(∇Nₐ, xₐ, ξ)
            ξ = [0.0, 0.0]
            ∇ₓNₐ = compute∇ₓNₐ(∇Nₐ(ξ), Jᵢⱼ(ξ))
            @test ∇ₓNₐ ≈ ∇Nₐ(ξ)

            # Scaling
            A = affineTransformationMatrix_2D(scale=[2,2])
            xₐ = buildLocalNodeCoordinates_2D(degree)
            for n = 1:size(xₐ,:local_node_id)
                xₐ[n,:] = (A * [xₐ[n,:]...,0.0])[1:2]
            end
            Jᵢⱼ = ξ->compute∇GeometricMapping(∇Nₐ, xₐ, ξ)
            ξ = [0.0, 0.0]
            ∇ₓNₐ = compute∇ₓNₐ(∇Nₐ(ξ), Jᵢⱼ(ξ))
            @test ∇ₓNₐ ≈ 1/2 * ∇Nₐ(ξ)
        end
        @testset "Degree = 2" begin
            degree = 2
            ∇Nₐ = ξ->∇LagrangeBasis_2D(degree,ξ)
            
            # Identity
            xₐ = buildLocalNodeCoordinates_2D(degree)
            Jᵢⱼ = ξ->compute∇GeometricMapping(∇Nₐ, xₐ, ξ)
            ξ = [0.0, 0.0]
            ∇ₓNₐ = compute∇ₓNₐ(∇Nₐ(ξ), Jᵢⱼ(ξ))
            @test ∇ₓNₐ ≈ ∇Nₐ(ξ)

            # Scaling
            A = affineTransformationMatrix_2D(scale=[2,2])
            xₐ = buildLocalNodeCoordinates_2D(degree)
            for n = 1:size(xₐ,:local_node_id)
                xₐ[n,:] = (A * [xₐ[n,:]...,0.0])[1:2]
            end
            Jᵢⱼ = ξ->compute∇GeometricMapping(∇Nₐ, xₐ, ξ)
            ξ = [0.0, 0.0]
            ∇ₓNₐ = compute∇ₓNₐ(∇Nₐ(ξ), Jᵢⱼ(ξ))
            @test ∇ₓNₐ ≈ 1/2 * ∇Nₐ(ξ)
        end
    end

    @testset "Dimension = 3" begin
        @testset "Degree = 1" begin
        end
        @testset "Degree = 2" begin
        end
    end
end

@testset "Newton-Raphson" begin
    @testset "1-D Linear" begin
        x̄ = newton_raphson(x->2-x, x->1, 10., 1, 1, 1e-12) 
        @test x̄ ≈ 2.0
    end

    @testset "1-D Quadratic" begin
        x̄ = newton_raphson(x->2-x^2, x->2*x, 10., 10, 10, 1e-12) 
        @test x̄ ≈ √2.0
    end

    @testset "1-D Cubic" begin
        x̄ = newton_raphson(x->2-x^3, x->3*x^2, 10., 10, 10, 1e-12) 
        @test x̄ ≈ (2.0)^(1/3)
    end
end

@testset "Nodal Coordinates" begin
    @testset "Dimension = 1" begin
        @testset "Degree = 1" begin
            deg = 1
            nodeCoords = [-1., 1.]
            @test buildLocalNodeCoordinates_1D(deg) ≈ nodeCoords
        end

        @testset "Degree = 2" begin
            deg = 2
            nodeCoords = [-1., 0., 1.]
            @test buildLocalNodeCoordinates_1D(deg) ≈ nodeCoords
        end
    end

    @testset "Dimension = 2" begin
        @testset "Degree = 1" begin
            deg = 1
            nodeCoords = [-1. -1.;
                           1. -1.;
                          -1.  1.;
                           1.  1.]
            @test buildLocalNodeCoordinates_2D(deg) ≈ nodeCoords
        end

        @testset "Degree = 2" begin
            deg = 2
            nodeCoords = [-1. -1.;
                           0. -1.;
                           1. -1.;
                          -1.  0.;
                           0.  0.;
                           1.  0.;
                          -1.  1.;
                           0.  1.;
                           1.  1.]
            @test buildLocalNodeCoordinates_2D(deg) ≈ nodeCoords
        end
    end
    
    @testset "Dimension = 3" begin
        @testset "Degree = 1" begin
            deg = 1
            nodeCoords = [-1. -1. -1.; 
                           1. -1. -1.; 
                          -1.  1. -1.;
                           1.  1. -1.; 
                          -1. -1.  1.;
                           1. -1.  1.; 
                          -1.  1.  1.;
                           1.  1.  1.]
            @test buildLocalNodeCoordinates_3D(deg) ≈ nodeCoords
        end

        @testset "Degree = 2" begin
            deg = 2
            nodeCoords = [-1. -1. -1.; 
                           0. -1. -1.;
                           1. -1. -1.;
                          -1.  0. -1.;
                           0.  0. -1.;
                           1.  0. -1.;
                          -1.  1. -1.;
                           0.  1. -1.;
                           1.  1. -1.;
                          -1. -1.  0.;
                           0. -1.  0.;
                           1. -1.  0.;
                          -1.  0.  0.;
                           0.  0.  0.;
                           1.  0.  0.;
                          -1.  1.  0.;
                           0.  1.  0.;
                           1.  1.  0.;
                          -1. -1.  1.;
                           0. -1.  1.;
                           1. -1.  1.;
                          -1.  0.  1.;
                           0.  0.  1.;
                           1.  0.  1.;
                          -1.  1.  1.;
                           0.  1.  1.;
                           1.  1.  1.]
            @test buildLocalNodeCoordinates_3D(deg) ≈ nodeCoords
        end
    end
end