using SparseArrays
using LinearAlgebra
using IncompleteLU
using Krylov
using LinearOperators

function constrPrecond(G⁻¹,A)
    #Fonction qui construit P⁻¹ pour un bloc G = diag(M)
    m, n = size(A)
    F2 = (A*G⁻¹*(A')) #Possibilite de faire des factorisation incomplete
    Aᵀ⁺GA⁺ = LinearOperator(Float64, m, m, true, true, u -> F2\u) 
    opM = LinearOperator(Float64, n+m, n+m, true, true, u -> [(G⁻¹-G⁻¹*A'*Aᵀ⁺GA⁺*A*G⁻¹)*u[1:n] + G⁻¹*A'*Aᵀ⁺GA⁺*u[n+1:m+n]; Aᵀ⁺GA⁺*A*G⁻¹*u[1:n]-Aᵀ⁺GA⁺*u[n+1:m+n]]) 
    return opM
end

A2 = sprand(100, 100, 0.9)
A1 = A2 + 0.00001*I
A1ILUT_1=IncompleteLU.ilu(A1,τ=1)
L = Matrix(A1ILUT_1.L+I)
U = Matrix(A1ILUT_1.U')
G⁻¹ = inv(L*U)
#println(norm(inv(Matrix(A1))-G⁻¹))
A = sprand(90,100,0.3)
N = spzeros(90,90)
m,n = size(A)
mat = [A1 A'; A -N]
mat2 = [A2 A'; A -N]
b=mat*ones(190)
opM = constrPrecond(G⁻¹,A)


res1 = zeros(3*m)
for k in 1:3*m
    x, stats = cgs(mat, b, itmax = k)
    res1[k] = norm(b - mat2*x)
end

res2 = zeros(3*m)
for k in 1:3*m
    x2, stats = cgs(mat, b, N = opM, itmax = k)
    res2[k] = norm(b - mat2*x2)
end

println("Residu du vecteur [x,y] du système par blocs sans préconditionneur: ", res1[3*m])
println("Residu du vecteur [x,y] du système par blocs avec préconditionneur: ", res2[3*m])

res = hcat(res1,res2)
using Plots
plot(1:3*m, res, yaxis=:log10, label=["sans préconditionneur" "avec préconditionneur iLU"])
plot!(ylabel = "residue ||rk||")
xlabel!("iteration k")
plot!(title="CGS")