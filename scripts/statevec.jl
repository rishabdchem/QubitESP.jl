using QubitESP 
const qesp = QubitESP 

#==================================
State vectors.
==================================#

let
    
    # Print
    println("Start the ESP state circuit")

    # Inputs 
    m = 8 
    n = 4 
    eta = randn(ComplexF64, m) 

    # ESP state
    U = espstate(eta, n) 
    V = focktohilbert(U, m, n) 

    # For printing 
    S1 = qesp.hamwtstrings(m, n)
    S2, H = qesp.qbitstrings(m)
    
    # Prepare output 
    f = open("results/output.txt", "w")

        println(f, "----------------------")

        println(f, "Number of qubits: $m, Hamming weight: $n")
        println(f, "Qubit indices are labelled from the right")

        println(f, "----------------------")

        println(f, )
        println(f, "ESP state vector in the qubit Fock space") 
        println(f, )

        for p = 1:2^m
            whatind = S2[p][length(S2[p])-m+1:end] 
            hamwt = H[p] 
            coeff = U[p] 
            println(f, "$whatind (Hamming wt: $hamwt): $coeff")
        end
        println(f, )

        println(f, "----------------------")

        println(f, )
        println(f, "ESP state vector in the qubit Hilbert space") 
        println(f, )

        for p = 1:binomial(m, n)
            whatind = S1[p]
            coeff = V[p] 
            println(f, "$whatind: $coeff")
        end
        println(f, )

        println(f, "----------------------")

        println(f, )
        println(f, "Scaled ESP state vector in the qubit Hilbert space") 
        println(f, )

        V = statescale!(V) 
        for p = 1:binomial(m, n)
            whatind = S1[p]
            coeff = V[p] 
            println(f, "$whatind: $coeff")
        end
        println(f, )

    close(f)

    # Print
    println("Check output file")

    # Exit    
    nothing    
end 
