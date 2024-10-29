using Test
using LinearAlgebra
include("../src/FewBodyPhysics.jl")

using .FewBodyPhysics

@testset "FewBodyPhysics.jl" begin

    @testset "Corput Sequence Tests" begin
        # Basic tests for base 2
        @test corput(1, 2) ≈ 0.5
        @test corput(2, 2) ≈ 0.25
        @test corput(3, 2) ≈ 0.75
        @test corput(4, 2) ≈ 0.125
        @test corput(5, 2) ≈ 0.625
    
        # Additional tests for base 3
        @test corput(1, 3) ≈ 1/3
        @test corput(2, 3) ≈ 2/3
        @test corput(3, 3) ≈ 1/9
        @test corput(4, 3) ≈ 4/9
        @test corput(5, 3) ≈ 7/9
    
        # Additional tests for base 5
        @test corput(1, 5) ≈ 0.2
        @test corput(2, 5) ≈ 0.04
        @test corput(3, 5) ≈ 0.6
        @test corput(4, 5) ≈ 0.08
        @test corput(5, 5) ≈ 0.16
        @test corput(10, 5) ≈ 0.32
    
    end
    
    @testset "Halton Sequence Tests" begin
        # Basic sequence tests for dimension 2
        @test halton(1, 2) ≈ [0.5, 1/3]
        @test halton(2, 2) ≈ [0.25, 2/3]
        @test halton(3, 2) ≈ [0.75, 1/9]
        @test halton(4, 2) ≈ [0.125, 4/9]
        @test halton(5, 2) ≈ [0.625, 7/9]
    
        # Additional tests for higher dimensions
        @test halton(1, 3) ≈ [0.5, 1/3, 0.2]
        @test halton(2, 3) ≈ [0.25, 2/3, 0.4]
        @test halton(3, 3) ≈ [0.75, 1/9, 0.6]
        @test halton(4, 3) ≈ [0.125, 4/9, 0.8]
        
        
        # Edge cases: checking very high `n`
        @test halton(1000, 2) ≈ [corput(1000, 2), corput(1000, 3)]
        @test halton(5000, 3) ≈ [corput(5000, 2), corput(5000, 3), corput(5000, 5)]
    
        # Edge case: requesting only 1 dimension
        @test halton(1, 1) ≈ [0.5]
        @test halton(10, 1) ≈ [corput(10, 2)]
    
        # Edge case: invalid dimension exceeding base list
        @test_throws AssertionError halton(1, 100)
    
        # Edge case: `n = 0` should return zero for any dimension
        @test halton(0, 2) ≈ [0.0, 0.0]
        @test halton(0, 5) ≈ [0.0, 0.0, 0.0, 0.0, 0.0]
    
        # Random large n tests for robustness
        @test halton(123, 4) ≈ [corput(123, 2), corput(123, 3), corput(123, 5), corput(123, 7)]
        @test halton(999, 3) ≈ [corput(999, 2), corput(999, 3), corput(999, 5)]
    end
    
    
end

