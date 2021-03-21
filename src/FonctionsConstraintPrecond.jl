using LinearAlgebra
using LinearOperators
using Krylov

## Fonctions qui créent différents types de G⁻¹ à partir de M
function blocGJacobi(M)

    #Preconditionneur avec bloc G diagonal   
    G⁻¹ = jacobi(M)
    return G⁻¹
end

## Fonctions qui "inversent" un certain G
function jacobi(A)

    # Inversion de G = diag(M)
    n = size(A, 1)
    J = zeros(n)
    for i = 1 : n
        J[i] = (A[i,i] == 0) ? 1.0 : 1 / A[i,i]
    end
    P⁻¹ = Diagonal(J)
    return P⁻¹
end

## Fonctions qui construit le préconditionneur
function constrPrecond(G⁻¹,A,N)

    #Calcul du complément de Schur et la méthodologie de son inversion
    m, n = size(A)
    F2 = Hermitian(N+A*G⁻¹*(A')) #Possibilite de faire des factorisations
    Schur_inv = LinearOperator(Float64, m, m, true, true, u -> F2\u) 

    opM = LinearOperator(Float64, n+m, n+m, true, true, u -> [(G⁻¹-G⁻¹*A'*Schur_inv*A*G⁻¹)*u[1:n] + G⁻¹*A'*Schur_inv*u[n+1:m+n]; Schur_inv*A*G⁻¹*u[1:n]-Schur_inv*u[n+1:m+n]]) 

    return opM

end

## Fonction qui résout le systeme global avec un préconditionneur par contrainte et une méthode donnée en entrée
function solvePrecond(M,A,N,methodKrylov="cgs", formG="Diagonal")

    #Construction de G⁻¹ du préconditionneur
    if formG == "Diagonal"
        G⁻¹ = blocGJacobi(M)
    else
        error("Cette forme de G n'est pas supportée")
    end

    #Construction du préconditionneur
    opM = constrPrecond(G⁻¹,A,N)

    mat = [M A'; A -N]

    #Différentes méthodes de Krylov possibles
    if methodKrylov == "cgs"
        x, stats = cgs(mat, b, itmax = 1000)
        x2, stats2 = cgs(mat, b, M = opM, itmax = 1000)
    elseif methodKrylov == "gmres"
        x, stats = dqgmres(mat, b, itmax = 1000)
        x2, stats2 = dqgmres(mat, b, M = opM, itmax = 1000)
    else
        error("Cette méthode de Krylov n'est pas supportée")
    end

    nbiter = length(stats.residuals) - 1
    println("Convergence en $nbiter itérations sans préconditionneur.")
    nbiter2 = length(stats2.residuals) - 1
    println("Convergence en $nbiter2 itérations sans préconditionneur.")

    println("Residu du vecteur [x,y] du système par blocs sans préconditionneur: ", norm(mat*x-b))
    println("Residu du vecteur [x,y] du système par blocs avec préconditionneur: ", norm(mat*x2-b))

    return x2, stats2

end