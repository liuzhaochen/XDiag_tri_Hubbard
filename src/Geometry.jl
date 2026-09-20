"""Rectangular quotient in triangular primitive-vector coordinates.

Coordinates are coefficients of `(a1,a2)`, not Cartesian coordinates. Sites are
ordered exactly as `site(x,y) = 1 + y + Ly*x` for zero-based coordinates.
"""
struct TriangularTorus
    Lx::Int
    Ly::Int
    function TriangularTorus(Lx::Integer, Ly::Integer)
        Lx >= 2 || throw(ArgumentError("Lx must be at least 2"))
        Ly >= 2 || throw(ArgumentError("Ly must be at least 2"))
        new(Int(Lx), Int(Ly))
    end
end

Base.length(torus::TriangularTorus) = torus.Lx * torus.Ly

@inline function site_index(torus::TriangularTorus, x::Integer, y::Integer)
    return 1 + mod(y, torus.Ly) + torus.Ly * mod(x, torus.Lx)
end

@inline function site_coordinate(torus::TriangularTorus, site::Integer)
    1 <= site <= length(torus) || throw(BoundsError(1:length(torus), site))
    offset = Int(site) - 1
    return (offset ÷ torus.Ly, offset % torus.Ly)
end

@inline _add_coordinate(a::NTuple{2,Int}, b::NTuple{2,Int}) =
    (a[1] + b[1], a[2] + b[2])

@inline _canonical_bond(i::Int, j::Int) = i < j ? (i, j) : (j, i)

"""All unique undirected nearest-neighbor bonds.

The positive displacements are `a1`, `a2`, and `a2-a1`.
"""
function nearest_neighbor_bonds(torus::TriangularTorus)
    positive_displacements = ((1, 0), (0, 1), (-1, 1))
    bonds = Set{NTuple{2,Int}}()
    for x in 0:(torus.Lx - 1), y in 0:(torus.Ly - 1)
        i = site_index(torus, x, y)
        for (dx, dy) in positive_displacements
            j = site_index(torus, x + dx, y + dy)
            i == j || push!(bonds, _canonical_bond(i, j))
        end
    end
    return sort!(collect(bonds))
end

function coordination_numbers(torus::TriangularTorus, bonds=nearest_neighbor_bonds(torus))
    coordination = zeros(Int, length(torus))
    for (i, j) in bonds
        coordination[i] += 1
        coordination[j] += 1
    end
    return coordination
end

function translation_permutation(torus::TriangularTorus, dx::Integer, dy::Integer)
    return [
        begin
            x, y = site_coordinate(torus, site)
            site_index(torus, x + dx, y + dy)
        end for site in 1:length(torus)
    ]
end

"""Squared Euclidean length of `x*a1+y*a2` in triangular-lattice units."""
@inline triangular_norm2(r::NTuple{2,Int}) = r[1]^2 + r[1] * r[2] + r[2]^2
