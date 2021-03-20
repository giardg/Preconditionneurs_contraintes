#using BlockArrays
using LinearAlgebra
using LinearOperators
using LimitedLDLFactorizations, IncompleteLU, ILUZero

function constrPrecond(M,A)

    #Fonction qui construit P⁻¹ pour un bloc G = diag(M)

    m, n = size(A)
    G = Diagonal(M)
    G⁻¹ = inv(G)
    F1 = A*G⁻¹*(A')
    #y = zeros(m,1);
    #Aᵀ⁺GA⁺ = LinearOperator(Float64, n, n, true, true, u -> F1\u) 

    opM = [I (-G⁻¹*(A')); zeros(m,n) I]*[G⁻¹ zeros(n,m); zeros(m,n) (-inv(F1))]*[I zeros(n,m);-A*G⁻¹ I] #a modifier bien sur pour ne pas avoir inv() et utilisé les LinearOperators

    #P = [G A'; A zeros(m,m)]

    return opM


end