using QubitESP
using Test
using Combinatorics 

@testset "ESP" begin

    # Initialize
    m = rand(4:20)
    n = rand(2:m-1)
    V = randn(m)
    println("m: $m, n: $n")

    # ESP
    s1 = espval(V, n) 

    # Brute force ESP
    X = Vector(1:m)
    comb = combinations(X, n)
    s2 = 0.0 
    for p in comb 
        t = 1.0 
        for j = 1:n
            t *= V[p[j]]  
        end
        s2 += t 
    end
    @test isapprox(s1, s2, atol=1e-12) 

end


@testset "On qubits" begin

    S1 = "10011"
    S2 = "100111"
    S3 = "1000110"
    @test QubitESP.findonbits(S1) == [1; 2; 5] 
    @test QubitESP.findonbits(S2) == [1; 2; 3; 6] 
    @test QubitESP.findonbits(S3) == [2; 3; 7] 

end


@testset "Fock to Hilbert" begin

    m = 2
    n = 1
    U = randn(2^m)
    V = focktohilbert(U, m, n)
    @test V == [U[2]; U[3]] 

end

