using ITensors

function make_index(names, values)
    return Tuple([Pair(name,value) for (name,value) in zip(names,values)])
end

function set_vector_helper(tensor, index, vector_index, values)
    tensor[(vector_index=>:, index...)...] = values
end

function pad_partition(partition, num_vars)
    new_partition = zeros(k)
    new_partition[1:length(partition)] = partition
    return new_partition
end

function pad_indices(permutation, num_active_vars, total_vars)
    new_index = zeros(total_vars)
    new_index[(total_vars - num_active_vars + 1):end] = new_index
    return new_index + ones(total_vars)
end

function generate_indices(monomial, num_active_vars, total_vars)
    if length(monomial) > num_active_vars
        return []
    padded_monomial = pad_partition(monomial, num_active_vars)
    permuted_monomials = multiset_permutations(padded_monomial, num_active_vars)
    return [pad_indices(index, num_active_vars, total_vars) for index in permuted_monomials]
end

# function generate_subterm(term)
#     nnz = count(x->x>0, term)
#     subterms = zeros(nnz, length(term))
#     index_counter = 0
#     for i=1:length(term)
#     end
# end