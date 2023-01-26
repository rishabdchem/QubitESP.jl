using LinearAlgebra 
using Combinatorics 

import LinearAlgebra.BLAS: BlasFloat

export espval, focktohilbert, statescale!

#====================================
ESP functions.
====================================#

"""
Compute ESP using the recursion:

S(m, n) = S(m-1, n) + V(m) S(m-1, n-1). 

# Arguments
- `V`: vector of dimension m 
- `n`: ESP degree 
"""
function espval(V::AbstractVector{T}, 
                n::Int) where T<:BlasFloat 

    # Check
    m = size(V, 1)
    m < 0 && error("m < 0 in espval")    
    n < 0 && error("n < 0 in espval")    

    # Trivial
    m < n && return zero(T) 
    n == 0 && return one(T) 
    n == 1 && return sum(V) 
    n == m && return prod(V) 

    # SumESP
    Smn = espval(V[1:m-1], n) + (V[m] * espval(V[1:m-1], n-1)) 

    return Smn 
end


"""
ESP factor: sqrt( S(p-1, q-1) / S(p, q) ).
Splitting is in the first index.  

# Arguments
- `V`: ESP vector of dimension p 
- `q`: degree index 
"""
function espfac(V::AbstractVector{T},
                q::Int) where T<:BlasFloat 

    # Check
    p = size(V, 1)
    p < q && error("p < q in espfac")
    p < 2 && error("p < 2 in espfac")
    q < 1 && error("q < 1 in espfac")

    # ESPs  
    f1 = espval(V[2:p], q-1)
    f2 = espval(V, q)

    return sqrt(f1 / f2) 
end

#======================================
Qubit state things.
======================================#

"""
Extract basis states with a given Hamming weight.

# Argument
- `U`: State vector 
- `m`: number of qubits 
- `n`: Hamming weight 
"""
function focktohilbert(U::AbstractVector{T},
                       m::Int,
                       n::Int) where T<:BlasFloat 

    # Check
    size(U, 1) == 2^m || error("Wrong input vector dimension in focktohilbert")

    # Initialize
    V = zeros(T, binomial(m, n))
    comb, M = mapvec(m, n)

    # Extract 
    a = 0 
    for p in comb 
        a += 1 
        V[a] = copy(U[M[a]])
    end

    return V
end


"""
<ESP|ESP>. 

# Arguments
- `eta`: ESP vector for m qubits
- `n`: Hamming weight 
"""
function espoverlap(eta::AbstractVector{T}, 
                    n::Int) where T<: BlasFloat
  
    X = zeros(T, size(eta, 1)) 
    X .= abs.(eta) .^ 2

    return espval(X, n) 
end


"""
Scale a state by multiplying with a global factor.
It does not change any observable.

# Arguments
- `V`: a state vector 
- `m`: number of qubits 
- `n`: Hamming weight 
"""
function statescale!(V::AbstractVector{T}) where T<: BlasFloat
 
    eps = 1e-04
    a = 0
    for j = 1:size(V, 1)
        a += 1 
        abs(V[j]) > eps && break   
        error("Check the input vector for threshold: $eps in statescale!")
    end 

    return V ./ V[a]
end

#======================================
Indexing stuffs.
======================================#

"""
Get one-bit positions from the right.

# Argument
- `s`: bitstring 
"""
function findonbits(s::String) 

    # Initialize
    indon = Int[]
    d = length(s)
    V = findall("1", s)
    n = size(V, 1)

    # Get locations
    for j = n:-1:1
        push!( indon, d - V[j][] + 1 ) 
    end  

    return indon 
end


"""
Ordering change from {0, 1}^m to {p1 < ... < pn}.

# Arguments
- `m`: number of qubits 
- `n`: Hamming weight 
"""
function mapvec(m::Int,
                n::Int)

    S, H = qbitstrings(m)
    V = zeros(Int, binomial(m, n))
    comb = combinations(Vector(1:m), n) 
    a = 0
    for p in comb
        a += 1
        for j = 1:2^m
            if p == findonbits(S[j])
                V[a] = j          
                break
            end 
        end
    end

    return comb, V
end


"""
Get qubit states as strings and their Hamming weights. 
Format: {0, 1}^m.

# Argument
- `m`: number of qubits 
"""
function qbitstrings(m::Int) 

    S = String[]
    H = Int[]
    for j = 1:2^m
        push!( S, bitstring(j - 1) )
        push!( H, count_ones(j - 1) ) 
    end

    return S, H
end


"""
Get qubit states as strings for a given Hamming weight n. 
Format: p1 < ... < pn. 

# Argument
- `m`: number of qubits 
- `n`: Hamming weight 
"""
function hamwtstrings(m::Int, 
                      n::Int) 

    S = String[]
    comb = combinations(Vector(1:m), n)
    for p in comb 
        str = ""
        for q = 1:m
            a = 0
            for j = 1:n
                if q == p[j] 
                    a = 1
                    str = "1" * str
                    break
                end 
            end 
            a == 1 && continue 
            str = "0" * str
        end
        push!( S, str) 
    end

    return S
end


