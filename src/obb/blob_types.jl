struct blob
    k::Int64
    n_below::Int64
    cost::Float64
    neighbor::Set{Int64}
    bin_BB_Tree::bin_BB_Tree
    function blob(k::Int64, n_below::Int64, neighbor_in, bb_tree::bin_BB_Tree, scale::Float64)
        bbb = bb_tree.box
        bbb = OBB(bbb, bbb)
        cost = blobCost(bbb, n_below, scale)
        return new(k, n_below, cost, Set{Int64}(neighbor_in), bb_tree)
    end
end

function createBlobDictionary(eM::eMesh, vec_obb_leaf::Vector{OBB}, vec_tri_tet, scale::Float64)
    vec_neighbor, is_abort = extractTriTetNeighborInformation(vec_tri_tet)
    dict_blob = Dict{Int64,blob}()
    if is_abort
        return dict_blob, true
    end
    for (k, ind_k) = enumerate(vec_tri_tet)
        aabb_k = vec_obb_leaf[k]
        tree_k = bin_BB_Tree(k, aabb_k)
        dict_blob[k] = blob(k, 1, vec_neighbor[k], tree_k, scale)
    end
    return dict_blob, false
end

function extractTriTetNeighborInformation(vec_tri_tet::Vector{SVector{N,Int64}}) where {N}
    dict_vert_pair, is_abort = createSharedEdgeFaceDict(vec_tri_tet)
    vec_neighbor = [Vector{Int64}() for _ = vec_tri_tet]
    if is_abort
        return vec_neighbor, true
    end
    for key_k = keys(dict_vert_pair)
        i_pair = dict_vert_pair[key_k]
        if (N == 3) || all(i_pair .!= -9999)  # tets will not have a partner if on the boundary
            push!(vec_neighbor[i_pair[1]], i_pair[2])
            push!(vec_neighbor[i_pair[2]], i_pair[1])
        end
    end
    return vec_neighbor, false
end

function createSharedEdgeFaceDict(vec_tri_tet::Vector{SVector{N,Int64}}) where {N}
    dict_vert_pair = Dict{SVector{N-1, Int64}, MVector{2,Int64}}()
    for (k_tet_tri, i_vertices) = enumerate(vec_tri_tet)
        for j = 1:N
            sv_new = sortEdgeFace(i_vertices, j)
            if haskey(dict_vert_pair, sv_new)
                mv_existing = dict_vert_pair[sv_new]
                if mv_existing[2] != -9999
                    (N == 3) && error("three triangles share the same edge something is wrong")
                    (N == 4) && error("three tetrahedrons share the same face something is wrong")
                end
                mv_existing[2] = k_tet_tri
            else
                dict_vert_pair[sv_new] = MVector{2,Int64}(k_tet_tri, -9999)
            end
        end
    end
    if N == 3  # confirm that every edge has a partner
        for (key_k, val_k) = dict_vert_pair
            if any(val_k .== -9999)
                #  && error("edge $key_k is connected to a single triangle $(val_k[1])")
                return dict_vert_pair, true
            end
        end
    end
    return dict_vert_pair, false
end

function blobCost(aabb::OBB, n_below::Int64, scale::Float64)
    cA = 1.0
    cV = 1.0
    V = 0.0
    V += n_below * log2(2 * n_below)
    V += cA * area(aabb) / (scale^2)
    V += cV * volume(aabb) / (scale^3)
    return V  # PriorityQueue returns lowest first
end

function doCombineBlob(dict_blob, a::blob, b::blob, k_next::Int64, scale::Float64)
    delete!(a.neighbor, b.k)
    delete!(b.neighbor, a.k)
    c_neighbor = union(a.neighbor, b.neighbor)  # add to neighbor_c
    tree_c = bin_BB_Tree(a.bin_BB_Tree, b.bin_BB_Tree)
    n_below_c = a.n_below + b.n_below
    blob_c = dict_blob[k_next] = blob(k_next, n_below_c, c_neighbor, tree_c, scale)
    return k_next + 1, blob_c
end

function calcMarginalCost(a::blob, b::blob, scale::Float64)
    c_cost = blobCost(OBB(a.bin_BB_Tree.box, b.bin_BB_Tree.box), a.n_below + b.n_below, scale)
    return c_cost - a.cost - b.cost
end

function refreshCostQueue!(dict_blob, pq_delta_cost, blob_a::blob, blob_c::blob, scale::Float64)
    a_k = blob_a.k
    for b_k = blob_a.neighbor  # for each neighbor in blob_a
        delete!(pq_delta_cost, minmax(b_k, a_k))  # delete old cost
        pq_delta_cost[minmax(b_k, blob_c.k)] = calcMarginalCost(dict_blob[b_k], blob_c, scale)  # add new cost
        blob_k = dict_blob[b_k]
        (a_k in blob_k.neighbor) || error("something is wrong")
        delete!(blob_k.neighbor, a_k)
        push!(blob_k.neighbor, blob_c.k)
    end
    delete!(dict_blob, a_k)
end

function createBlobPriorityQueue(dict_blob, scale)
    pq_delta_cost = PriorityQueue{Tuple{Int64,Int64},Float64}()
    for k_a = keys(dict_blob)
        blob_a = dict_blob[k_a]
        for k_b = blob_a.neighbor
            (k_a < k_b) && (pq_delta_cost[(k_a, k_b)] = calcMarginalCost(blob_a, dict_blob[k_b], scale))
        end
    end
    return pq_delta_cost
end

function bottomUp!(dict_blob, pq_delta_cost, scale)
    k_next = length(dict_blob) + 1
    while !isempty(pq_delta_cost)
        key_lowest, val_lowest = dequeue_pair!(pq_delta_cost)  # automatically removes entry
        a_k, b_k = key_lowest
        blob_a = dict_blob[a_k]
        blob_b = dict_blob[b_k]
        k_next, blob_c = doCombineBlob(dict_blob, blob_a, blob_b, k_next, scale)
        refreshCostQueue!(dict_blob, pq_delta_cost, blob_a, blob_c, scale)  # deal with blob_a
        refreshCostQueue!(dict_blob, pq_delta_cost, blob_b, blob_c, scale)  # deal with blob_b
    end
    return nothing
end

function eMesh_to_tree(eM::eMesh{T1,T2}) where {T1,T2}
    isa(eM, eMesh{Tri,Tet}) && error("Cannot create tree for eMesh{Tri,Tet} use as_tri_eMesh or as_tet_eMesh on input first.")
    n_leaf = ifelse(T1==Tri, n_tri, n_tet)(eM)
    if (n_leaf == 1)
        if T1 == Tri
            obb = calc_obb(eM.point[eM.tri[1]])
        else
            obb = calc_obb(eM.point[eM.tet[1]])
        end
        return bin_BB_Tree(1, obb)
    end
    point = eM.point
    vec_tri_tet = ifelse(isa(eM, eMesh{Tri,Nothing}), eM.tri, eM.tet)
    scale = sum(calc_obb(point).e) / 3
    if T1 == Tri
        vec_obb_leaf = [calc_obb(eM.point[k]) for k = eM.tri]
    else
        vec_obb_leaf = [calc_obb(eM.point[k]) for k = eM.tet]
    end
    dict_blob, is_abort = createBlobDictionary(eM, vec_obb_leaf, vec_tri_tet, scale)
    is_abort && error("not implemented error: disconnected mesh")

    pq_delta_cost = createBlobPriorityQueue(dict_blob, scale)
    bottomUp!(dict_blob, pq_delta_cost, scale)
    all_tree = collect(values(dict_blob))
    if (length(all_tree) != 1)
        error("not implemented error: disconnected mesh")
    else
        tree = all_tree[1].bin_BB_Tree
        tight_fit_leaves!(eM, tree)
        return tree
    end
end

function tight_fit_leaf!(eM::eMesh{Tri,Nothing}, eT::bin_BB_Tree{OBB})
    eT.box = fit_tri_obb(vertex_pos_for_tri_ind(eM, eT.id))
end

function tight_fit_leaf!(eM::eMesh{Nothing,Tet}, eT::bin_BB_Tree{OBB})
    eT.box = fit_tet_obb(vertex_pos_for_tet_ind(eM, eT.id), eM.ϵ[eM.tet[eT.id]])
end

function tight_fit_leaves!(eM::eMesh, eT::bin_BB_Tree{OBB})
    if is_leaf(eT)
        tight_fit_leaf!(eM, eT)
    else
        tight_fit_leaves!(eM, eT.node_1)
        tight_fit_leaves!(eM, eT.node_2)
    end
end
