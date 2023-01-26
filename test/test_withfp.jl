#=========================================================
This is an extra test related to the connection
between AGP and ESP state.
It needs the FermiPairing.jl package:
github.com/rishabdchem/FermiPairing.jl. 
=========================================================#

using QubitESP
using Test
using FermiPairing #github.com/rishabdchem/FermiPairing.jl 
const fp = FermiPairing

@testset "AGP energy for pairing Hamiltonian" begin

    # Initialize
    m = rand(4:12)
    n = rand(2:m-1)
    eta = rand(ComplexF64, m)
    G = randn()
    println("m: $m, n: $n")

    # Test
    V = focktohilbert( espstate(eta, n), m, n ) 
    e1 = docihamoverlap(V, V, m, n, G) / docioverlap(V, V, m, n)
    e2 = agphamoverlap(eta, eta, n, G) / agpoverlap(eta, eta, n) 
    @test isapprox(e1, e2, atol=1e-08) 

end

