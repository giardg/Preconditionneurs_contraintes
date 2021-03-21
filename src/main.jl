using SparseArrays
include("FonctionsConstraintPrecond.jl")

A = rand(90,120)
M = rand(120,120)
N = sprand(90,90,0.1)
M = M'*M #M = M'*M pour hermitien, mais pas necessaire
N = N'*N;

m,n = size(A)
b = rand(m+n)

#Resolution avec un preconditionneur par contrainte (G = diag(M) par défaut)
x, stats = solvePrecond(M,A,N);
