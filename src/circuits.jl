using Yao 

import LinearAlgebra.BLAS: BlasFloat

#==================================
ESP circuits.
==================================#

"""
ESP magnitude circuit.

# Argument
- `theta`: rotation matrix 
"""
function maguni(theta::AbstractMatrix{T}) where T<:BlasFloat 

    # Initialize
    m = size(theta, 1) 
    n = size(theta, 2)

    # Build
    C = scsgate(theta[m, 1:n], (m, n), m, 1)
    a = 1
    for p = m-1:-1:n+1
        a += 1 
        push!( C, scsgate(theta[p, 1:n], (p, n), m, a) )
    end
    for p = n:-1:2
        a += 1 
        push!( C, scsgate(theta[p, 1:p-1], (p, p-1), m, a) )
    end
   
    return C 
end


"""
Split and shift cyclic circuit. 

# Arguments
- `theta`: rotation vector 
- `gind`: gate indices
- `m`: circuit width 
- `qind`: first qubit index 
"""
function scsgate(theta::AbstractVector{T},
                 gind::Tuple,
                 m::Int,
                 qind::Int) where T<:BlasFloat 

    # Check
    size(theta, 1) == gind[2] || error("Dimension mismatch in scsgate")

    # Initialize
    p, q = gind[1], gind[2]

    # Build
    C = twogate(2acos(theta[1]), m, (qind, qind+1,))
    for r = 2:q         
        push!( C, threegate(2acos(theta[r]), m, (qind, qind+r-1, qind+r,)) )
    end
   
    return C 
end


"""
Two-qubit gate.

# Arguments
- `tau`: Ry angle 
- `m`: circuit width 
- `qind`: qubit indices
"""
function twogate(tau::T,
                 m::Int,
                 qind::Tuple) where T<:BlasFloat 

    # Check
    size(qind, 1) == 2 || error("Wrong number of indices in twogate")
    qind[1] > m && error("Wrong index in twogate")
    qind[2] > m && error("Wrong index in twogate")

    # Build
    p, q = qind[1], qind[2]
    G = chain(m,
        control(q, p=>X),
        control(p, q=>Ry(tau)),
        control(q, p=>X))

    return G 
end


"""
Three-qubit gate.

# Arguments
- `tau`: Ry angle 
- `m`: circuit width 
- `qind`: qubit indices
"""
function threegate(tau::T,
                   m::Int,
                   qind::Tuple) where T<:BlasFloat 

    # Check
    size(qind, 1) == 3 || error("Wrong number of indices in threegate")
    qind[1] > m && error("Wrong index in threegate")
    qind[2] > m && error("Wrong index in threegate")
    qind[3] > m && error("Wrong index in threegate")

    # Build
    p, q, r = qind[1], qind[2], qind[3]
    G = chain(m,
        control(r, p=>X),
        control((p, q,), r=>Ry(tau)),
        control(r, p=>X))

    return G 
end


