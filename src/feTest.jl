using Test
import Base.≈
≈(a::Float64,b::Float64) = isapprox(a::Float64,b::Float64;atol=10*eps(),rtol=10*eps())

include("feCode.jl")
include("feElement.jl")

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
        nodeCoords = [[-1., -1.], [1., -1.], [-1., 1.],[1., 1.]]
        for n = 1:numNodes
            @test LagrangeBasis_2D(degree, nodeCoords[n])[n] ≈ 1.0
        end
    end

    @testset "Degree = 2" begin
        degree = 2
        numNodes = (degree + 1)^2
        nodeCoords = [[-1., -1.], [0., -1.], [1., -1.], [-1., 0.],[0., 0.], [1., 0.], [-1., 1.], [0., 1.] ,[1., 1.]]
        for n = 1:numNodes
            @test LagrangeBasis_2D(degree, nodeCoords[n])[n] ≈ 1.0
        end
    end
end

# Test 3D Lagrange Basis Function
@testset "3D Lagrange Basis Function" begin
    @testset "Degree = 1" begin
        degree = 1
        numNodes = (degree + 1) ^ 3
        nodeCoords = [[-1., -1., -1.], [1., -1., -1.], 
                      [-1.,  1., -1.], [1.,  1., -1.], 
                      [-1., -1.,  1.], [1., -1.,  1.], 
                      [-1.,  1.,  1.], [1.,  1.,  1.]]
        for n = 1:numNodes
            @test LagrangeBasis_3D(degree, nodeCoords[n])[n] ≈ 1.0
        end
    end

    @testset "Degree = 2" begin
        degree = 2
        numNodes = (degree + 1)^3
        nodeCoords = [[-1., -1., -1.], [0., -1., -1.], [1., -1., -1.], 
                      [-1.,  0., -1.], [0.,  0., -1.], [1.,  0., -1.],
                      [-1.,  1., -1.], [0.,  1., -1.], [1.,  1., -1.],
                      [-1., -1.,  0.], [0., -1.,  0.], [1., -1.,  0.], 
                      [-1.,  0.,  0.], [0.,  0.,  0.], [1.,  0.,  0.],
                      [-1.,  1.,  0.], [0.,  1.,  0.], [1.,  1.,  0.],
                      [-1., -1.,  1.], [0., -1.,  1.], [1., -1.,  1.], 
                      [-1.,  0.,  1.], [0.,  0.,  1.], [1.,  0.,  1.],
                      [-1.,  1.,  1.], [0.,  1.,  1.], [1.,  1.,  1.]]
        for n = 1:numNodes
            @test LagrangeBasis_3D(degree, nodeCoords[n])[n] ≈ 1.0
        end
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