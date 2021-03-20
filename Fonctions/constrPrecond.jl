#using BlockArrays
using LinearAlgebra
using LinearOperators
using LimitedLDLFactorizations, IncompleteLU, ILUZero

function blocGJacobi(M)

    #Preconditionneur avec bloc G diagonal   
    G⁻¹ = jacobi(M)
    return G⁻¹
end

function constrPrecond(G⁻¹,A)

    #Fonction qui construit P⁻¹ pour un bloc G = diag(M)

    m, n = size(A)
    F2 = (A*G⁻¹*(A')) #Possibilite de faire des factorisation incomplete
    Aᵀ⁺GA⁺ = LinearOperator(Float64, m, m, true, true, u -> F2\u) 

    opM = LinearOperator(Float64, n+m, n+m, true, true, u -> [(G⁻¹-G⁻¹*A'*Aᵀ⁺GA⁺*A*G⁻¹)*u[1:n] + G⁻¹*A'*Aᵀ⁺GA⁺*u[n+1:m+n]; Aᵀ⁺GA⁺*A*G⁻¹*u[1:n]-Aᵀ⁺GA⁺*u[n+1:m+n]]) 

    return opM


end

function jacobi(A)
    n = size(A, 1)
    J = zeros(n)
    for i = 1 : n
        J[i] = (A[i,i] == 0) ? 1.0 : 1 / A[i,i]
    end
    P⁻¹ = Diagonal(J)
    return P⁻¹
end