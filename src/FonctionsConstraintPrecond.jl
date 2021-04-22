using LinearOperators
using Krylov
using LimitedLDLFactorizations,LDLFactorizations, IncompleteLU, ILUZero

## Fonctions qui créent différents types de G⁻¹ à partir de M (Fonctions inutilisées)
function blocGJacobi(M::AbstractArray)

    #Preconditionneur avec bloc G diagonal   
    G⁻¹ = jacobi(M)
    return G⁻¹
end

function lim_LDL(A)
    n = size(A, 1)
    F = lldl(A, memory=10)
    F.D = abs.(F.D);
    G⁻¹ =  LinearOperator(Float64, n, n, true, true, v -> F \ v) 
end

## Fonctions qui "inversent" un certain G (Fonction inutilisée)
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

## Fonction qui construit le préconditionneur à gauche ou à droite (Fonction inutilisée)
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
function solvePrecond(M::AbstractArray,A::AbstractArray,N::AbstractArray,D,methodKrylov="gmres", precond = true, formG="Diagonal", precondGDC="G", rtol=0.0, atol=0.0, maxit = 1000)

    m,n = size(A)
    
    if precond
        #Construction de G⁻¹ et de G du préconditionneur
        if formG == "Diagonal"
            G = Diagonal(M)

        elseif formG == "Symmetric"
            G = (1/2).*(M+M') #Partie symétrique de M

        elseif formG == "I"
            G = LinearOperator(Float64, n, n, true, true, v -> v) 
        else
            error("Cette forme de G n'est pas supportée")
        end

        #Construction du préconditionneur (par inversion par bloc) 
        #opM = constrPrecond(G⁻¹,A,N) #Probablement <a retirer

        #Construction  du préconditionneur (par factorisation lldl)
        P = ([Matrix(G) A'; A -N])

        if precondGDC == "G"
            F = lldl(P,memory=1000)
            opM = LinearOperator(Float64, m+n, m+n, true, true, w -> F\w)
            opN = LinearOperator(Float64, m+n, m+n, true, true, v -> v) 
        elseif precondGDC == "D" 
            F = lldl(P,memory=1000)
            opN = LinearOperator(Float64, m+n, m+n, true, true, w -> F\w)
            opM = LinearOperator(Float64, m+n, m+n, true, true, v -> v) 
        elseif precondGDC=="C"
            F = lu(P)
            opM = LinearOperator(Float64, m+n, m+n, false, false, u -> F.L\u)
            opN = LinearOperator(Float64, m+n, m+n, false, false, v -> F.U\v)
        else
            error("Le préconditionneur ne peut être appliqué qu'à gauche (G), droite (D) ou centré (C)")
        end

    else
        #P = I si on ne veut pas de preconditionneur
        opM = LinearOperator(Float64, m+n, m+n, true, true, v -> v) 
        opN = LinearOperator(Float64, m+n, m+n, true, true, v -> v) 
    end

    mat = [M A'; A -N]

    #Différentes méthodes de Krylov possibles
    if methodKrylov == "cgs"
        x, stats = cgs(mat, D, M = opM, N = opN, rtol=rtol, atol=atol, itmax=maxit,history=true)
    elseif methodKrylov == "gmres"
        x, stats = dqgmres(mat, D, M = opM, N = opN, rtol=atol, atol=atol, itmax=maxit,history=true)
    elseif methodKrylov == "bicgstab"
        x, stats = bicgstab(mat, D, M = opM, N = opN, rtol=rtol, atol=atol, itmax=maxit,history=true)
    elseif methodKrylov == "diom"
        x, stats = diom(mat, D, M = opM, N = opN, rtol=rtol, atol=atol, itmax=maxit,memory=1000,history=true)
    else
        error("Cette méthode de Krylov n'est pas supportée")
    end

    return x, stats

end
