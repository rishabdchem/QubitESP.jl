using Yao 
using Combinatorics 

import LinearAlgebra.BLAS: BlasFloat

export espstate, dickestate

#=======================================================================
We want to prepare the ESP state

|Psi> = sum(1 <= p1 < ... < pn <= m) eta(p1) ... eta(pn) 
        | ... 0 ... 1(pn) ... 1(p1) ... 0 ... >,

where all eta(p) = |eta(p)| exp(i phi(p)). 
When all eta(p) = 1, it is called the Dicke state.

We have two formats for qubit basis state ordering

Fock space: {0, 1}^m,
Hilbert space: 1 <= p1 < ... < pn <= m,

where m = number of qubits and n = Hamming weight.
=======================================================================#

"""
Build ESP state in the Fock space.

# Arguments
- `eta`: ESP state of m qubits 
- `n`: Hamming weight 
"""
function espstate(eta::AbstractVector{T},
                  n::Int) where T<:BlasFloat 

    # Check
    n >= size(eta, 1) && error("n >= m in espstate") 
    
    # Initialize
    m = size(eta, 1) 
    norm = zeros(T, m)
    theta = zeros(T, m, n)
    etamag = zeros(Float64, m)
    phi = zeros(Float64, m)
    etamag .= abs.(Complex.(eta))
    phi .= angle.(Complex.(eta))
 
    # Rotation angles  
    norm .= etamag .^ 2
    for q = 1:n, p = q+1:m 
        a = m - p + 1
        theta[p, q] = etamag[a] * espfac(norm[a:m], q)  
    end

    # Initial circuit
    C1 = chain(m, repeat(X, 1:n))  

    # Magnitude circuit 
    C2 = maguni(theta)

    # Phase circuit 
    C3 = chain(m, put(p=>Rz(phi[p])) for p = 1:m)
    
    # Final state
    psi = zero_state(m) |> C1 |> C2 |> C3  

    return statevec(psi)
end


"""
Build Dicke state in the Fock space.

# Arguments
- `m`: number of qubits 
- `n`: Hamming weight 
"""
function dickestate(m::Int,
                    n::Int) 

    # Check
    n >= m && error("n >= m in dickestate") 
    
    # Initialize
    theta = zeros(Float64, m, n)
 
    # Rotation angles  
    for q = 1:n, p = q+1:m 
        a = m - p + 1
        theta[p, q] = sqrt(q / p) 
    end

    # Initial circuit
    C1 = chain(m, repeat(X, 1:n))  

    # Magnitude circuit 
    C2 = maguni(theta)

    # Final state
    psi = zero_state(m) |> C1 |> C2 

    return statevec(psi)
end


"""
Build ESP state in the Hilbert space with brute force. 

# Arguments
- `eta`: ESP state of m qubits 
- `n`: Hamming weight 
"""
function bfespstate(eta::AbstractVector{T},
                    n::Int) where T<:BlasFloat 

    # Check
    n > size(eta, 1) && error("n > m in bfespstate") 

    # Initialize
    m = size(eta, 1) 
    d = binomial(m, n)
    V = ones(T, d)

    # Get combinations
    X = Vector(1:m)
    comb = combinations(X, n)

    # DOCI coefficients 
    mu = 0
    for p in comb
        mu += 1
        for j = 1:n    
            V[mu] *= eta[p[j]]
        end    
    end 

    return Complex.(V) 
end

