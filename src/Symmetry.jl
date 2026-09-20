@inline function _matrix_on_coordinate(A::Matrix{Int}, r::NTuple{2,Int})
    return (A[1, 1] * r[1] + A[1, 2] * r[2],
            A[2, 1] * r[1] + A[2, 2] * r[2])
end

@inline _matrix_key(A::Matrix{Int}) = (A[1, 1], A[1, 2], A[2, 1], A[2, 2])

"""Generate the 12 integer-coordinate operations of triangular-lattice `D6`."""
function integer_point_group()
    operations = Matrix{Int}[]
    for a in -2:2, b in -2:2, c in -2:2, d in -2:2
        abs(a * d - b * c) == 1 || continue
        col1, col2 = (a, c), (b, d)
        triangular_norm2(col1) == 1 || continue
        triangular_norm2(col2) == 1 || continue
        triangular_norm2(_add_coordinate(col1, col2)) == 3 || continue
        push!(operations, [a b; c d])
    end
    sort!(operations; by=_matrix_key)
    length(operations) == 12 || error("failed to generate triangular D6")
    return operations
end

function _preserves_torus(torus::TriangularTorus, A::Matrix{Int})
    # M=diag(Lx,Ly); the quotient is preserved iff M^-1*A*M is integral.
    return mod(A[1, 2] * torus.Ly, torus.Lx) == 0 &&
           mod(A[2, 1] * torus.Lx, torus.Ly) == 0
end

"""Subgroup of triangular `D6` that preserves the finite torus."""
function surviving_point_group(torus::TriangularTorus)
    return [A for A in integer_point_group() if _preserves_torus(torus, A)]
end

function point_permutation(torus::TriangularTorus, A::Matrix{Int})
    return [
        begin
            r = site_coordinate(torus, site)
            x, y = _matrix_on_coordinate(A, r)
            site_index(torus, x, y)
        end for site in 1:length(torus)
    ]
end

function preserves_bonds(
    torus::TriangularTorus,
    A::Matrix{Int},
    bonds=nearest_neighbor_bonds(torus),
)
    permutation = point_permutation(torus, A)
    transformed = Set(
        _canonical_bond(permutation[i], permutation[j]) for (i, j) in bonds
    )
    return transformed == Set(bonds)
end

"""Check identity, inverses, and multiplication closure for site permutations."""
function is_closed_permutation_group(permutations::Vector{Vector{Int}})
    isempty(permutations) && return false
    nsites = length(first(permutations))
    keys = Set(Tuple(permutation) for permutation in permutations)
    Tuple(1:nsites) in keys || return false
    all(sort(permutation) == collect(1:nsites) for permutation in permutations) ||
        return false
    for left in permutations, right in permutations
        Tuple(left[right[index]] for index in 1:nsites) in keys || return false
    end
    return true
end
