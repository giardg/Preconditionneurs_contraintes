using Krylov
using SparseArrays
include("constrPrecond.jl")

A = rand(90,120)
M = rand(120,120)
M = M'*M

m,n = size(A)

mat = [M A'; A zeros(m,m)]
b = rand(m+n)

G = blocMPrecond(M)
opM = constrPrecond(G,A)

x, stats = dqgmres(mat, b)
nbiter = length(stats.residuals) - 1
println("Convergence en $nbiter itérations sans préconditionneur.")

x2, stats = dqgmres(mat, b, M = opM)
nbiter = length(stats.residuals) - 1
println("Convergence en $nbiter itérations avec préconditionneur.")

println("Residu du vecteur [x,y] du système par blocs sans préconditionneur: ", norm(mat*x-b))
println("Residu du vecteur [x,y] du système par blocs avec préconditionneur: ", norm(mat*x2-b))
