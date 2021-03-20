using Krylov
include("constrPrecond.jl")

A = rand(90,120)
M = rand(120,120)
M = M'*M

m,n = size(A)

mat = [M A'; A zeros(m,m)]
b = rand(m+n)

opM = constrPrecond(M,A)

println(typeof(mat))
x, stats = dqgmres(mat, b)
nbiter = length(stats.residuals) - 1
println("Convergence en $nbiter itérations sans préconditionneur.")

mat2 = opM*mat
b2 = opM*b
x2, stats = dqgmres(mat2, b2)
nbiter = length(stats.residuals) - 1
println("Convergence en $nbiter itérations avec préconditionneur.")

println("Residu du vecteur [x,y] du système par blocs sans préconditionneur: ", norm(mat*x-b))
println("Residu du vecteur [x,y] du système par blocs avec préconditionneur: ", norm(mat*x2-b))
