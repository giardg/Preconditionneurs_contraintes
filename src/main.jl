using SparseArrays
using Plots
include("FonctionsConstraintPrecond.jl")

#Test sur des systemes randoms pour commencer
A = sprand(90,120,0.1)
N = sprand(90,90,0.0) #N = 0 souvent, mais pas nécessaire
M = sprand(120,120,0.1)
M = M'*M #M = M'*M pour hermitien, mais pas nécessaire (pour gmres oui)
N = N'*N;

m,n = size(A)
b = [rand(m);rand(n)]
mat = [M A'; A -N]

#Resolution avec un preconditionneur par contrainte (G = diag(M) par défaut)
res1 = zeros(100)
res2 = zeros(100)
for k = 1:100
    x, stats, x2, stats2 = solvePrecond(M,A,N,b, "gmres", "Diagonal", k);
    res1[k] = norm(b-mat*x)
    res2[k] = norm(b-mat*x2)
end

plot(1:100, res1, yaxis=:log10, xlabel="k", color=4, ylabel="|r_k|", label="Sans préconditionneur")
plot!(1:100, res2, yaxis=:log10, xlabel="k", color=5, ylabel="|r_k|", label="Avec préconditionneur")
