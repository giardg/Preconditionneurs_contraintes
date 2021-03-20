# Aim of this code :

# The goal of this code is to know is the matrix N intervening in the saddle point matrix
# [M A'; A -N] is semi-definite positive or not.

# Theorical explanation :

# All Hermitian matrices are diagonalizable (Spectral theorem). An Hermitian matrix is positive 
# semi-definite iff its eigenvalues are greater or equal to 0. So, in order to know if Hermitian 
# matrix is positive semi-definite, it is sufficient to look at whether its eigenvalues are 
# positive or zero thanks to Julia's function eigvals().

# Libraries :

using SparseArrays, LinearAlgebra

# 1st case : A symmetric matrix (dense).

function test_definite_semi_positive(M::Array{Float64,2})
    count = 0
    eigenvalues_M = eigvals(M)
    if typeof(eigenvalues_M) == Array{Complex{Float64},1}
        return false
    end
    for k in 1:length(eigenvalues_M)
        if eigenvalues_M[k]>=0
            count+=1
        end
    end
    if count == length(eigenvalues_M)
        return true
    else 
        return false
    end
end

# 1st case bis : A symmetric matrix (sparse).

function test_definite_semi_positive(M::SparseMatrixCSC{Float64,Int64})
    count = 0
    eigenvalues_M = eigen(collect(M)).values
    if typeof(eigenvalues_M) == Array{Complex{Float64},1}
        return false
    end
    for k in 1:length(eigenvalues_M)
        if eigenvalues_M[k]>=0
            count+=1
        end
    end
    if count == length(eigenvalues_M)
        return true
    else 
        return false
    end
end


# 2nd case : A hermitian matrix (dense).

function test_definite_semi_positive(M::Array{Complex{Float64},2})
    count = 0
    eigenvalues_M = eigvals(M)
    if typeof(eigenvalues_M) == Array{Complex{Float64},1}
        return false
    end
    for k in 1:length(eigenvalues_M)
        if eigenvalues_M[k]>=0
            count+=1
        end
    end
    if count == length(eigenvalues_M)
        return true
    else 
        return false
    end
end

# 2nd case bis : A hermitian matrix (sparse).

function test_definite_semi_positive(M::SparseMatrixCSC{Complex{Float64},Int64})
    count = 0
    eigenvalues_M = eigen(collect(M)).values
    if typeof(eigenvalues_M) == Array{Complex{Float64},1}
        return false
    end
    for k in 1:length(eigenvalues_M)
        if eigenvalues_M[k]>=0
            count+=1
        end
    end
    if count == length(eigenvalues_M)
        return true
    else 
        return false
    end
end