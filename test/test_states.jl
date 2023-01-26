using QubitESP
using Test

@testset "Dicke state" begin

    # Initialize
    m = rand(4:12)
    n = rand(2:m-1)
    println("m: $m, n: $n")

    # Test
    U = focktohilbert( dickestate(m, n), m, n )
    V = ones(ComplexF64, binomial(m, n)) 
    @test isapprox(statescale!(U), statescale!(V), atol=1e-12) 

end


@testset "ESP state" begin

    # Initialize
    m = 4 #rand(4:12)
    n = 2 #rand(2:m-1)
    eta = randn(ComplexF64, m)
    println("m: $m, n: $n")

    # Test
    U = focktohilbert( espstate(eta, n), m, n )
    V = QubitESP.bfespstate(eta, n)
    @test isapprox(statescale!(U), statescale!(U), atol=1e-12) 

end


@testset "Check intruder qubit states" begin

    # Initialize
    m = rand(4:12)
    n = rand(2:m-1)
    eta = randn(ComplexF64, m)
    println("m: $m, n: $n")

    # ESP state
    V = espstate(eta, n)
    S, H = QubitESP.qbitstrings(m)

    # Test
    eps = 1e-12
    a = 0
    for j = 1:2^m
        count_ones(j - 1) == n && continue
        if abs(V[j]) > eps
            a = 1 
            break
        end   
    end
    @test a == 0 

end


