function pad_partition(partition, total_vars)
    new_partition = zeros(total_vars)
    new_partition[1:length(partition)] = partition
    return new_partition
end

function pad_indices(permutation, num_active_vars, total_vars)
    new_index = zeros(total_vars)
    new_index[(total_vars - num_active_vars + 1):end] = permutation
    return Tuple(new_index)
end

function generate_indices(monomial, num_active_vars, total_vars)
    if length(monomial) > num_active_vars
        return []
    end
    padded_monomial = pad_partition(monomial, num_active_vars)
    permuted_monomials = multiset_permutations(padded_monomial, num_active_vars)
    return [pad_indices(index, num_active_vars, total_vars) for index in permuted_monomials]
end

function increment_index(index, var)
    update = collect(index)
    update[var] += 1
    return Tuple(update)
end

# function generate_subterm(term)
#     nnz = count(x->x>0, term)
#     subterms = zeros(nnz, length(term))
#     index_counter = 0
#     for i=1:length(term)
#     end
# end