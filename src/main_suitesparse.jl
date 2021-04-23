using SparseArrays, MAT, LinearAlgebra
using SuiteSparseMatrixCollection, MatrixMarket
using Plots
using LaTeXStrings
include("FonctionsConstraintPrecond.jl")

matrix_name = "garon1"
pb = ssmc_matrices("", matrix_name)
fetch_ssmc(pb, format="MM")
path_mtx = matrix_paths(pb, format="MM")[1]
mat = MatrixMarket.mmread(joinpath(path_mtx, matrix_name * ".mtx"));

n=size(mat,1)
for i = 1:size(mat,1)
   if mat[i,i] == 0
        global n = i #On assume qu'un problème avec N = 0 a été choisi
        break
   end
end


#n = 14 #Pour bfwb62 qui n'a pas N=0 

m = size(mat,1)-n;
M = mat[1:n,1:n]
A = mat[n+1:m+n,1:n]
N = -mat[n+1:m+n,n+1:m+n]
D = [rand(n); rand(m)]

println(norm(N'-N))


t1 = @elapsed x1, stats = solvePrecond(M,A,N,D, "gmres", false, "I", "G", 1e-8, 1e-8, 50000);
println(matrix_name, " -- ", " Sans préconditionneur Nombre d'itérations: ", length(stats.residuals)-1, " Temps: ", t1)
println(matrix_name, " -- ", " Sans préconditionneur Résidu: ", norm((mat*x1-D)))

t2 = @elapsed x2, stats = solvePrecond(M,A,N,D, "gmres", true, "I", "G", 1e-8, 1e-8, 50000);
println(matrix_name, " -- ", " G = I Nombre d'itérations: ", length(stats.residuals)-1, " Temps: ", t2)
println(matrix_name, " -- ", " G = I Résidu: ", norm(mat*x2-D))

t3 = @elapsed x3, stats = solvePrecond(M,A,N,D, "gmres", true, "Diagonal", "G", 1e-8, 1e-8, 50000);
println(matrix_name, " -- ", " G = Diag(M) Nombre d'itérations: ", length(stats.residuals)-1, " Temps: ", t3)
println(matrix_name, " -- ", " G = Diag(M) Résidu: ", norm(mat*x3-D))

t4 = @elapsed x4, stats = solvePrecond(M,A,N,D, "gmres", true, "Symmetric", "G", 1e-8, 1e-8, 50000);
println(matrix_name, " -- ", " G = 1/2(M+M') Nombre d'itérations: ", length(stats.residuals)-1, " Temps: ", t4)
println(matrix_name, " -- ", " G = 1/2(M+M') Résidu: ", norm(mat*x4-D))

kmax = 100
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

#savefig(plot_ref1,"C:\\Users\\Gregory.Giard\\Maitrise\\Session 2\\Algebre lineaire numerique appliquee\\Projet\\Preconditionneurs_contraintes\\gmres_gauche_rk_Maxwell")