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