using LinearAlgebra
using LinearOperators
using Krylov
using LimitedLDLFactorizations

## Fonctions qui créent différents types de G⁻¹ à partir de M
function blocGJacobi(M::AbstractArray)

    #Preconditionneur avec bloc G diagonal   
    G⁻¹ = jacobi(M)
    return G⁻¹
end

function lim_LDL(A)
    n = size(A, 1)
    F = lldl(A, memory=1)
    F.D = abs.(F.D);
    G⁻¹ =  LinearOperator(Float64, n, n, true, true, v -> F \ v) 
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

## Fonction qui construit le préconditionneur à gauche ou à droite
function constrPrecond(G⁻¹,A::AbstractArray,N::AbstractArray)

    #Calcul du complément de Schur et la méthodologie de son inversion
    m, n = size(A)
    F2 = Hermitian(Matrix(N+A*G⁻¹*(A'))) #Possibilite de faire des factorisations ici
    F2 = cholesky(F2);
    Schur_inv = LinearOperator(Float64, m, m, true, true, u -> F2\u) 

    opM = LinearOperator(Float64, n+m, n+m, true, true, u -> [(G⁻¹-G⁻¹*A'*Schur_inv*A*G⁻¹)*u[1:n] + G⁻¹*A'*Schur_inv*u[n+1:m+n]; Schur_inv*A*G⁻¹*u[1:n]-Schur_inv*u[n+1:m+n]]) 

    return opM

end

## Fonction qui résout le systeme global avec un préconditionneur par contrainte et une méthode donnée en entrée
function solvePrecond(M::AbstractArray,A::AbstractArray,N::AbstractArray,b,methodKrylov="cgs", formG="Diagonal", maxit = 1000)

    #Construction de G⁻¹ du préconditionneur
    if formG == "Diagonal"
        G⁻¹ = blocGJacobi(M)

    elseif formG == "LLDL"
        G⁻¹ = lim_LDL(M)

    elseif formG == "I"
        G⁻¹ = LinearOperator(Float64, size(M,1), size(M,1), true, true, v -> v) 
    else
        error("Cette forme de G n'est pas supportée")
    end

    #Construction du préconditionneur
    opM = constrPrecond(G⁻¹,A,N)

    mat = [M A'; A -N]

    #Différentes méthodes de Krylov possibles
    if methodKrylov == "cgs"
        x, stats = cgs(mat, b, rtol=0.0, atol=0.0, itmax=maxit,history=true)
        x2, stats2 = cgs(mat, b, M = opM, rtol=0.0, atol=0.0, itmax=maxit,history=true)
    elseif methodKrylov == "gmres"
        x, stats = dqgmres(mat, b, rtol=0.0, atol=0.0, itmax=maxit,history=true)
        x2, stats2 = dqgmres(mat, b, M = opM, rtol=0.0, atol=0.0, itmax=maxit,history=true)
    else
        error("Cette méthode de Krylov n'est pas supportée")
    end

    #nbiter = length(stats.residuals) - 1
    #println("Convergence en $nbiter itérations sans préconditionneur.")
    #nbiter2 = length(stats2.residuals) - 1
    #println("Convergence en $nbiter2 itérations sans préconditionneur.")

    #println("Residu du vecteur [x,y] du système par blocs sans préconditionneur: ", norm(mat*x-b))
    #println("Residu du vecteur [x,y] du système par blocs avec préconditionneur: ", norm(mat*x2-b))

    return x, stats, x2, stats2

end
