using SparseArrays, MAT, LinearAlgebra
using SuiteSparseMatrixCollection, MatrixMarket
include("FonctionsConstraintPrecond.jl")

matrix_name = "spaceShuttleEntry_4"
pb = ssmc_matrices("", matrix_name)
fetch_ssmc(pb, format="MM")
path_mtx = matrix_paths(pb, format="MM")[1]
mat = MatrixMarket.mmread(joinpath(path_mtx, matrix_name * ".mtx"));
println(size(mat))

n=size(mat,1)
for i = 1:size(mat,1)
   if mat[i,i] == 0
        global n = i
        break
   end
end

m = size(mat,1)-n;
M = mat[1:n,1:n]
A = mat[n+1:m+n,1:n]
N = -mat[n+1:m+n,n+1:m+n]
D = [rand(n); rand(m)]

x1, stats = solvePrecond(M,A,N,D, "gmres", false, "I", 1e-8, 1e-8, m+n);
println(matrix_name, " -- ", " Sans préconditionneur Nombre d'itérations: ", length(stats.residuals)-1)
println(matrix_name, " -- ", " Sans préconditionneur Résidu: ", norm((mat*x1-D)))

x2, stats = solvePrecond(M,A,N,D, "gmres", true, "I", 1e-8, 1e-8, m+n);
println(matrix_name, " -- ", " G = I Nombre d'itérations: ", length(stats.residuals)-1)
println(matrix_name, " -- ", " G = I Résidu: ", norm(mat*x2-D))

x3, stats = solvePrecond(M,A,N,D, "gmres", true, "Diagonal", 1e-8, 1e-8, m+n);
println(matrix_name, " -- ", " G = Diag(M) Nombre d'itérations: ", length(stats.residuals)-1)
println(matrix_name, " -- ", " G = Diag(M) Résidu: ", norm(mat*x3-D))
