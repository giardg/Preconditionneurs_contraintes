using SparseArrays, MAT, LinearAlgebra
using Plots
using LaTeXStrings
include("FonctionsConstraintPrecond.jl")

## Construction du système à résoudre

#Test sur des systemes randoms pour commencer
A = sprand(90,120,1.0)
N = sprand(90,90,0.0) #N = 0 souvent, mais pas nécessaire
M = sprand(120,120,1.0)
M = M'*M #M = M'*M pour hermitien, mais pas nécessaire (pour gmres oui)
N = N'*N;

#A[end,:] = A[1,:]; #Matrice A pas de rang plein

m,n = size(A)
mat = [M A'; A -N]
#D = [rand(n);zeros(m)] 
D = [rand(n);rand(m)]
#D = mat[:,1]


## Resolution avec un preconditionneur par contrainte
kmax = 500
res1 = zeros(kmax)
res2 = zeros(kmax)
res3 = zeros(kmax)
res4 = zeros(kmax)
time1 = zeros(kmax)
time2 = zeros(kmax)
time3 = zeros(kmax)
time4 = zeros(kmax)

for k = 1:kmax
    t1 = @elapsed x, stats = solvePrecond(M,A,N,D, "gmres", false, "I", "G", 0.0, 0.0, k);
    res1[k] = norm(D-mat*x)
    time1[k] = t1
    t2 = @elapsed x2, stats = solvePrecond(M,A,N,D, "gmres", true, "I", "G", 0.0, 0.0, k);
    res2[k] = norm(D-mat*x2)
    time2[k] = t2
    t3 = @elapsed x3, stats = solvePrecond(M,A,N,D, "gmres", true, "Diagonal", "G", 0.0, 0.0, k);
    res3[k] = norm(D-mat*x3)
    time3[k] = t3
    t4 = @elapsed x4, stats = solvePrecond(M,A,N,D, "gmres", true, "Symmetric", "G", 0.0, 0.0, k);
    res4[k] = norm(D-mat*x4)
    time4[k] = t4
end

## Affichage des résultats
plot_ref1=plot(1:kmax, res1, yaxis=:log10, xlabel=L"k", color=1, ylabel=L"||r_k||", label="Sans préconditionneur",legend=:best, fmt = :png)
plot_ref1=plot!(1:kmax, res2, color=2, label=L"G = I")
plot_ref1=plot!(1:kmax, res3, color=3, label=L"G = Diag(M)")
plot_ref1=plot!(1:kmax, res4, color=4, label=L"G = \frac{1}{2}(M+M^*)")

#savefig(plot_ref1,"C:\\Users\\Gregory.Giard\\Maitrise\\Session 2\\Algebre lineaire numerique appliquee\\Projet\\Preconditionneurs_contraintes\\gmres_gauche_rk")

