# Explanation :
# All Hermitian matrices are diagonalizable (Spectral theorem). An Hermitian matrix is positive 
# semi-definite iff its eigenvalues are greater or equal to 0. So, in order to know if Hermitian 
# matrix is positive semi-definite, it is sufficient to look at whether its eigenvalues are 
# positive or zero thanks to Julia's function eigvals().

# 1st case : A symmetric matrix.

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

# 2nd case : A hermitian matrix.

function test_definite_positive(M::Array{Complex{Float64},2})
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