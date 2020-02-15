using Test
import Base.≈
≈(a::Float64,b::Float64) = isapprox(a::Float64,b::Float64;atol=10*eps(),rtol=10*eps())

include("feCode.jl")

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