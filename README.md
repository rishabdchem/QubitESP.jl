# QubitESP.jl

State preparation of elementary symmetric polynomial (ESP) state. 
Algorithm implementation using the [Yao.jl](https://github.com/QuantumBFS/Yao.jl) framework.

## Introduction 

ESP state is a generalization of Dicke state, where instead of all coefficients being the same, they have an ESP structure.
Using fermion-qubit mapping, it can be shown that ESP is equivalent to a fermionic state known as antisymmetrized geminal power (AGP) in chemistry and number-projected Bardeen-Cooper-Schrieffer (PBCS) in physics.

## References
   
1. Armin Khamoshi, Rishab Dutta, Gustavo E. Scuseria. 
   State preparation of AGP on a quantum computer without number projection
   ([preprint](https://arxiv.org/abs/2301.09586)).
 
1. Andreas Bartschi, Stephen Eidenbenz. 
   Deterministic preparation of Dicke states.
   *Fundamentals of Computation Theory.* **2019**, Springer, Cham 
   ([preprint](https://arxiv.org/abs/1904.07358))
   ([article](https://doi.org/10.1007/978-3-030-25027-0_9)).
 
