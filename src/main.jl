using SparseArrays, MAT, LinearAlgebra
using Plots
using LaTeXStrings
include("FonctionsConstraintPrecond.jl")

## Construction du système à résoudre

#Test sur des systemes randoms pour commencer
A = sprand(90,120,1.0)
N = sprand(90,90,0.5) #N = 0 souvent, mais pas nécessaire
M = sprand(120,120,1.0)
M = M'*M #M = M'*M pour hermitien, mais pas nécessaire (pour gmres oui)
N = N'*N;

#A[end,:] = A[1,:]; #Matrice A pas de rang plein

m,n = size(A)
mat = [M A'; A -N]
D = [rand(m);zeros(n)] 
#D = [rand(m);rand(n)]
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
    #t4 = @elapsed x4, stats = solvePrecond(M,A,N,b, "gmres", true, "LLDL", k);
    #res4[k] = norm(b-mat*x4)
    #time4[k] = t4
end

## Affichage des résultats
plot_ref1=plot(1:kmax, res1, yaxis=:log10, xlabel=L"k", color=1, ylabel=L"\|r_k\|", label="Sans préconditionneur",legend=:best, fmt = :png)
plot_ref1=plot!(1:kmax, res2, color=2, label="G = I")
plot_ref1=plot!(1:kmax, res3, color=3, label="G = Diagonal(M)")
#plot_ref1=plot!(1:kmax, res4, yaxis=:log10, xlabel=L"k", color=4, ylabel=L"\|r_k\|", label="G = LLDL",legend=:best, fmt = :png)

#savefig(plot_ref1,"C:\\Users\\Gregory.Giard\\Maitrise\\Session 2\\Algebre lineaire numerique appliquee\\Projet\\Preconditionneurs_contraintes\\cgs_rk_nonHerm")

#plot_ref2=plot(1:kmax, time1, yaxis=:log10, xlabel=L"k", color=1, ylabel="Temps (s)", label="Sans préconditionneur",legend=:bottomright, fmt = :png)
#plot_ref2=plot!(1:kmax, time2, color=2, label="G = I")
#plot_ref2=plot!(1:kmax, time3, color=3, label="G = Diagonal(M)")
#plot_ref2=plot!(1:kmax, time4, yaxis=:log10, xlabel=L"k", color=4, ylabel="Temps (s)", label="G = LLDL",legend=:bottomright, fmt = :png)

#savefig(plot_ref2,"C:\\Users\\Gregory.Giard\\Maitrise\\Session 2\\Algebre lineaire numerique appliquee\\Projet\\Preconditionneurs_contraintes\\cgs_time_nonHerm")