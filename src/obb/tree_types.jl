mutable struct bin_BB_Tree{BB_Type}
    id::Int64
    box::BB_Type
    node_1::bin_BB_Tree
    node_2::bin_BB_Tree
    function bin_BB_Tree(id::Int64, bb::BB_Type) where {BB_Type <: BoundingBox}
        return new{BB_Type}(id, bb)
    end
    function bin_BB_Tree(node_1::bin_BB_Tree{BB_Type}, node_2::bin_BB_Tree{BB_Type}) where {BB_Type <: BoundingBox}
        aabb = BB_Type(node_1.box, node_2.box)
        return new{BB_Type}(-9999, aabb, node_1, node_2)
    end
    function bin_BB_Tree(bb::BB_Type, node_1::bin_BB_Tree{BB_Type}, node_2::bin_BB_Tree{BB_Type}) where {BB_Type <: BoundingBox}
        return new{BB_Type}(-9999, bb, node_1, node_2)
    end
end

function Base.show(io::IO, bin_tree::bin_BB_Tree)
    n_leaf = leafNumber(bin_tree)
    n_depth = treeDepth(bin_tree)
    println(io, "bin_BB_Tree with $n_leaf leaves in $n_depth layers:")
    if n_depth == 1
        frac_high = n_leaf / (2^n_depth - 1)
        frac_low = n_leaf / (2^(n_depth - 1) - 1)
        println("    fill between $(frac_high) and $(frac_low) of optimal")
        if (frac_high < 0.1) || (frac_low < 0.2)
            @warn("tree may be unbalanced")
        end
    end
end

mutable struct TT_Cache
    vc::VectorCache{Tuple{Int64, Int64}}
    t_a_b::SVector{3,Float64}
    R_a_b::SMatrix{3,3,Float64,9}
    abs_R_a_b::SMatrix{3,3,Float64,9}
    function TT_Cache()
        vc = VectorCache{NTuple{2,Int64}}()
        return new(vc)
    end
end

update_TT_Cache!(tt::TT_Cache, trans::SVector{3,Float64}, R::RotMatrix{3,Float64,9}) = update_TT_Cache!(tt, trans, SMatrix{3,3,Float64,9}(R))
function update_TT_Cache!(tt::TT_Cache, trans::SVector{3,Float64}, R::SMatrix{3,3,Float64,9})
    tt.t_a_b = trans
    tt.R_a_b = R
    tt.abs_R_a_b = abs.(R) .+ 1.0e-14
    empty!(tt.vc)
    return nothing
end

Base.length(tt::TT_Cache) = length(tt.vc)
Base.isempty(tt::TT_Cache) = isempty(tt.vc)
Base.@propagate_inbounds Base.getindex(tt::TT_Cache, i::Int) = tt.vc[i]

is_leaf(tree::bin_BB_Tree)     = (tree.id != -9999)
is_not_leaf(tree::bin_BB_Tree) = (tree.id == -9999)

treeDepth(t::bin_BB_Tree) = leafNumberDepth(t, 0, 0)[1]
leafNumber(t::bin_BB_Tree) = leafNumberDepth(t, 0, 0)[2]

function leafNumberDepth(t::bin_BB_Tree, k_depth::Int64, k_leaf::Int64)
    if is_leaf(t)
        k_leaf = 1
    else is_not_leaf(t)
        k_depth_1, k_leaf_1 = leafNumberDepth(t.node_1, k_depth + 1, k_leaf)
        k_depth_2, k_leaf_2 = leafNumberDepth(t.node_2, k_depth + 1, k_leaf)
        k_depth = max(k_depth_1, k_depth_2)
        k_leaf = k_leaf_1 + k_leaf_2
    end
    return k_depth, k_leaf
end

function extractData!(tree::bin_BB_Tree, v::Vector{Int64})
    if is_leaf(tree)
        push!(v, tree.id)
    else
        extractData!(tree.node_1, v)
        extractData!(tree.node_2, v)
    end
end
function extractData(tree::bin_BB_Tree)
    v = Vector{Int64}()
    extractData!(tree, v)
    return v
end

function tree_tree_intersect(ttCache::TT_Cache, tree_1::bin_BB_Tree, tree_2::bin_BB_Tree)
    if BB_BB_intersect(ttCache, tree_1.box, tree_2.box)
        is_leaf_1 = is_leaf(tree_1)
        is_leaf_2 = is_leaf(tree_2)
        if is_leaf_1
            if is_leaf_2
                addCacheItem!(ttCache.vc, (tree_1.id, tree_2.id))
            else
                tree_tree_intersect(ttCache, tree_1, tree_2.node_1)
                tree_tree_intersect(ttCache, tree_1, tree_2.node_2)
            end
        else
            if is_leaf_2
                tree_tree_intersect(ttCache, tree_1.node_1, tree_2)
                tree_tree_intersect(ttCache, tree_1.node_2, tree_2)
            else
                tree_tree_intersect(ttCache, tree_1.node_1, tree_2.node_1)
                tree_tree_intersect(ttCache, tree_1.node_2, tree_2.node_1)
                tree_tree_intersect(ttCache, tree_1.node_1, tree_2.node_2)
                tree_tree_intersect(ttCache, tree_1.node_2, tree_2.node_2)
            end
        end
    end
end
