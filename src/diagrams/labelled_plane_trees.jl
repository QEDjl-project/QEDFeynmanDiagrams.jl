function _gen_trees(vec::Vector{Tuple{Int64,Vector{Int64}}})
    # input is one tree structure
    # node neighbours are partially filled from the front

    n = length(vec)

    # recursion ends if only two nodes are left with 1 
    sum_open_edges = sum(getindex.(vec, 1)) # total number of open edges
    sum_open_edges -= sum(length.(getindex.(vec, 2))) # minus total number of assigned edges

    @assert sum_open_edges != 0 "cannot have empty tree"
    @assert sum_open_edges % 2 == 0 "cannot have uneven number of open edges"

    # find first entry that has only one open edge
    single_open_edge_nodes = Vector{Int64}()
    i = 0
    for node in vec
        i += 1
        if node[1] - length(node[2]) == 1
            push!(single_open_edge_nodes, i)
        end
    end

    @assert length(single_open_edge_nodes) >= 2 "there were less than two nodes with only one open edge, which cannot happen for valid trees"

    results = Vector{Vector{Vector{Int64}}}()

    if sum_open_edges == 2
        # assign open edges
        n1 = single_open_edge_nodes[1]
        n2 = single_open_edge_nodes[2]
        push!(vec[n1][2], n2)
        push!(vec[n2][2], n1)

        push!(results, getindex.(vec, 2))

        return results
    end

    # choose first node
    n1 = single_open_edge_nodes[1]

    # iterate through second nodes: all nodes that have at least 
    for n2 in 1:n
        if vec[n2][1] - length(vec[n2][2]) <= 1
            # must have at least 2 open slots to be a partner for the edge
            continue
        end

        # make new vec and connect edges
        recurse_vec = deepcopy(vec)
        push!(recurse_vec[n1][2], n2)
        push!(recurse_vec[n2][2], n1)

        generated_trees = _gen_trees(recurse_vec)
        append!(results, generated_trees)
    end

    return results
end

function _labelled_plane_trees_unpermuted(N::Int64)
    all_trees = Vector{Vector{Vector{Int64}}}()

    parts = partitions(N - 2)
    if N == 1
        push!(all_trees, Vector{Vector{Vector{Int64}}}())
        return all_trees
    end

    if N == 2
        # fix empty partition when n = 2
        parts = [Vector{Int64}()]
    end

    for origp in parts
        p = copy(origp) # can't change original partition because of the iteration

        n = length(p)
        resize!(p, N)
        for i in 1:N
            if i > n
                p[i] = 0
            end
            p[i] += 1
        end

        for mp in multiset_permutations(p, N)
            # representation as vector of vectors
            # each vector is a node
            # each vector has the indices of the connected nodes
            tree = [(l, Vector{Int64}()) for l in mp]

            generated_trees = _gen_trees(tree)
            append!(all_trees, generated_trees)
        end
    end

    return all_trees
end

function _plane_tree_permutations_helper(t::Vector{Vector{Int64}}, depth::Int64)
    if depth > length(t)
        return [t]
    end

    all_perms = Vector{Vector{Vector{Int64}}}()
    for permutation in permutations(t[depth])
        permuted_t = deepcopy(t)
        permuted_t[depth] = permutation
        append!(all_perms, _plane_tree_permutations_helper(permuted_t, depth + 1))
    end

    return all_perms
end

# return all the possible edge permutations of a given labelled plane tree
function _plane_tree_permutations(t::Vector{Vector{Int64}})
    return _plane_tree_permutations_helper(t, 1)
end

function labelled_plane_trees(N::Int64)
    if N == 1
        return [[Int[]]]
    end

    all_trees = Vector{Vector{Vector{Int64}}}()

    sizehint!(all_trees, factorial(3N - 3, 2N - 1))

    for t in _labelled_plane_trees_unpermuted(N)
        append!(all_trees, _plane_tree_permutations(t))
    end

    return all_trees
end
