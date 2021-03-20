using Krylov
using SparseArrays
include("constrPrecond.jl")

A = sprand(90,120,0.5)
M = sprand(120,120,0.5)
N = sprand(90,90,0.0)
M = M'*M #M = M'*M pour hermitien, mais pas necessaire
N = N'*N;

m,n = size(A)

mat = [M A'; A -N]
b = rand(m+n)

G⁻¹ = blocGJacobi(M)
opM = constrPrecond(G⁻¹,A,N)

x, stats = cgs(mat, b, itmax = 1000)
nbiter = length(stats.residuals) - 1
println("Convergence en $nbiter itérations sans préconditionneur.")

x2, stats = cgs(mat, b, N = opM, itmax = 1000)
nbiter = length(stats.residuals) - 1
println("Convergence en $nbiter itérations avec préconditionneur.")

println("Residu du vecteur [x,y] du système par blocs sans préconditionneur: ", norm(mat*x-b))
println("Residu du vecteur [x,y] du système par blocs avec préconditionneur: ", norm(mat*x2-b))
