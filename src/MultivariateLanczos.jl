using LinearAlgebra, Combinatorics
include("utils.jl")

"""
    MultivariateLanczos(A::Vector{Matrix{Float64}}, b::Vector{Float64}, n::Int64, m::Int64)

    Inputs-
        A : An ordered list of n x n support matrices Ai. Requires all matrices A1,...,Ak be Hermitian and simultaneously diagonalizable (i.e. commuting).

        b : Starting weight vector.

        n : Dimension of support matrices Ai.

        m : Maximum polynomial degree desired.

    Outputs-
        T : Dict{Tuple,Float64} indexed by ((i1,i2,...,ik), (j1,j2,...,jk)) such that T[((i1,i2,...,ik), (j1,j2,...,jk))] equals the recurrence coefficent in the multivariate orthogonal polynomial basis for the polynomial with leading term s1^(i1) * s2^(i2) * ... * sk^(ik) in terms of the polynomial with leading term s1^(j1) * s2^(j2) * ... * sk^(jk). The variables s1,...,sk are ordered by the input list A.

        Q : Dict{Tuple,Vector} indexed by (i1,i2,...,ik) such that T[(i1, ..., ik)] equals the orthogonal polynomial with leading term s1^(i1) * s2^(i2) * ... * sk^(ik) evaluated over the support (spectrum of matrices Ai).
"""
function MultivariateLanczos(A,b,n,m)
    TOL = 0.0001 # tolerance for orthogonality
    k = length(A) # number of variables

    T = Dict{Tuple,Float64}()
    Q = Dict{Tuple,Vector}()

    zero_index = Tuple(zeros(Int32, k))
    T[zero_index] = norm(b)
    Q[zero_index] = b / norm(b)

    for i=1:m # i = total degree
        for var=1:k
            num_active_vars = k - var + 1
            prev_monomials = i > 1 ? integer_partitions(i-1) : [[0]]
            prev_indices = [generate_indices(monom, num_active_vars,k) for monom in prev_monomials]

            for index_set in prev_indices # index set contains all monomials with same partition
                for index in index_set
                    if !haskey(Q,index)
                        continue
                    end

                    v = Q[index]
                    q = A[var] * v
                    q = q/norm(q)
                    next_index = increment_index(index, var)

                    # orthogonalize q against all previous vectors, can be optimized
                    for (vec_index, vec) in pairs(Q)
                        t = (vec' * q)
                        T[(vec_index, next_index)] = t
                        q = q - t * vec
                    end
                    if norm(q) > TOL
                        T[next_index] = norm(q)
                        Q[next_index] = q / norm(q)
                    end
                end
            end
        end
    end
    return (T,Q)
end
