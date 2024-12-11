using LinearAlgebra, ITensors, Combinatorics
include("utils.jl")

"""
    MultivariateLanczos(A::Dict{String, Matrix{Float64}}, b::Vector{Float64}, n::Int64, m::Int64)

    Inputs-
        A : {Pair(s1, A1), Pair(s2, A2), ..., Pair(sk, Ak)} is a dictionary pairing user-defined names for each variable si and its corresponding n x n support matrix Ai. Requires all matrices A1,...,Ak be Hermitian and simultaneously diagonalizable (i.e. commuting). Also requires that no variable name equals "MOP_value".

        b : Starting weight vector.

        n : Dimension of support matrices Ai.

        m : Maximum polynomial degree desired.

    Outputs-
        I : Dictionary mapping variable names to ITensor indices 

        T : ITensor indexed by (s1,s2,...,sk) such that T(s1=>i1, ..., sk=>ik) equals the recurrence coefficent in the multivariate orthogonal polynomial basis for the polynomial with leading term s1^(i1-1) * s2^(i2-1) * ... * sk^(ik-1).

        Q : ITensor indexed by (s1,s2,...,sk, MOP_value) such that T(s1=>i1, ..., sk=>ik, MOP_value => 1:n) equals the orthogonal polynomial with leading term s1^(i1-1) * s2^(i2-1) * ... * sk^(ik-1) evaluated over the support (spectrum of matrices Ai).
"""
function MultivariateLanczos(A,b,n,m)
    VECTOR_DIM_NAME = "MOP_value"

    k = length(A) # number of variables
    names_ordered = keys(A)
    tensor_indices = Tuple([Index(m+1,i) for i in names_ordered])
    vector_index = Index(n, VECTOR_DIM_NAME)
    indices_with_vector = (vector_index, tensor_indices...)

    T = ITensor(tensor_indices...)
    Q = ITensor(indices_with_vector...)

    zero_index = make_index(tensor_indices, ones(Int32, k))
    Q[(vector_index=>1, zero_index...)...] = 0.0 # need tensor to be nonempty
    T[zero_index...] = norm(b)
    Q[(vector_index=>:, zero_index...)...] = b / norm(b)

    for i=1:m # i = total degree
        for var=1:k
            num_active_vars = k - var + 1
            prev_monomials = i > 0 ? integer_partitions(i-1) : [[]]
            prev_indices = [generate_indices for degs in prev_monomials]

            for index_set in prev_indices # index set contains all monomials with same partition
                for index in index_set
                    tensor_index = make_index(tensor_indices, index)
                    v = 

                    T[zero_index...] = norm(b)
                    Q[(vector_index=>:, zero_index...)...] = b / norm(b)
                end
            end
        end
    end
    return Q
end
