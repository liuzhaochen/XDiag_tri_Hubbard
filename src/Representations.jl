struct MomentumSector
    nx::Int
    ny::Int
end

Base.:(==)(a::MomentumSector, b::MomentumSector) = a.nx == b.nx && a.ny == b.ny
Base.hash(k::MomentumSector, h::UInt) = hash((k.nx, k.ny), h)

function all_momenta(torus::TriangularTorus)
    momenta = MomentumSector[]
    for nx in 0:(torus.Lx - 1), ny in 0:(torus.Ly - 1)
        push!(momenta, MomentumSector(nx, ny))
    end
    return momenta
end

@inline _fractional_part(x::Rational{Int}) = x - floor(Int, x)

function _momentum_coordinates(torus::TriangularTorus, momentum::MomentumSector)
    return (momentum.nx // torus.Lx, momentum.ny // torus.Ly)
end

function _momentum_from_coordinates(
    torus::TriangularTorus,
    k::NTuple{2,Rational{Int}},
)
    k1, k2 = _fractional_part(k[1]), _fractional_part(k[2])
    nx = k1 * torus.Lx
    ny = k2 * torus.Ly
    denominator(nx) == 1 && denominator(ny) == 1 ||
        error("point operation did not preserve the allowed momentum grid")
    return MomentumSector(mod(numerator(nx), torus.Lx), mod(numerator(ny), torus.Ly))
end

"""Active reciprocal-space action `k -> A^(-T) k`, modulo reciprocal vectors."""
function _momentum_action(
    torus::TriangularTorus,
    A::Matrix{Int},
    momentum::MomentumSector,
)
    a, b, c, d = A[1, 1], A[1, 2], A[2, 1], A[2, 2]
    determinant = a * d - b * c
    abs(determinant) == 1 || throw(ArgumentError("point operation must be unimodular"))
    k1, k2 = _momentum_coordinates(torus, momentum)
    transformed = ((d * k1 - c * k2) / determinant,
                   (-b * k1 + a * k2) / determinant)
    return _momentum_from_coordinates(torus, transformed)
end

function momentum_star(torus::TriangularTorus, momentum::MomentumSector)
    return sort!(
        collect(Set(_momentum_action(torus, A, momentum)
                    for A in surviving_point_group(torus)));
        by=k -> (k.nx, k.ny),
    )
end

function production_momentum_representatives(torus::TriangularTorus)
    seen = Set{MomentumSector}()
    representatives = MomentumSector[]
    for momentum in all_momenta(torus)
        momentum in seen && continue
        push!(representatives, momentum)
        union!(seen, momentum_star(torus, momentum))
    end
    return representatives
end

struct AffineOperation
    matrix::Matrix{Int}
    displacement::NTuple{2,Int}
    permutation::Vector{Int}
end

function _affine_permutation(
    torus::TriangularTorus,
    A::Matrix{Int},
    displacement::NTuple{2,Int},
)
    return [
        begin
            r = _matrix_on_coordinate(A, site_coordinate(torus, site))
            site_index(torus, r[1] + displacement[1], r[2] + displacement[2])
        end for site in 1:length(torus)
    ]
end

function _translation_phase(
    torus::TriangularTorus,
    momentum::MomentumSector,
    displacement::NTuple{2,Int},
)
    phase = momentum.nx * displacement[1] / torus.Lx +
            momentum.ny * displacement[2] / torus.Ly
    return cis(2pi * phase)
end

function translation_operations_and_characters(
    torus::TriangularTorus,
    momentum::MomentumSector,
)
    identity_matrix = [1 0; 0 1]
    operations = AffineOperation[]
    characters = ComplexF64[]
    for x in 0:(torus.Lx - 1), y in 0:(torus.Ly - 1)
        displacement = (x, y)
        push!(operations, AffineOperation(
            identity_matrix,
            displacement,
            translation_permutation(torus, x, y),
        ))
        push!(characters, _translation_phase(torus, momentum, displacement))
    end
    return operations, characters
end

function _little_point_group(torus::TriangularTorus, momentum::MomentumSector)
    return [A for A in surviving_point_group(torus)
            if _momentum_action(torus, A, momentum) == momentum]
end

"""Little-space-group representation, with optional `C2` parity.

For this cluster a parity is valid only at `(0,0)` and `(0,2)`. Generic momenta
use translation symmetry alone in production.
"""
function little_group_operations_and_characters(
    torus::TriangularTorus,
    momentum::MomentumSector,
    parity::Int,
)
    parity in (-1, 1) || throw(ArgumentError("C2 parity must be +1 or -1"))
    little_group = _little_point_group(torus, momentum)
    length(little_group) == 2 ||
        throw(ArgumentError("momentum does not have the C2 little group"))
    identity_matrix = [1 0; 0 1]
    operations = AffineOperation[]
    characters = ComplexF64[]
    for A in little_group, x in 0:(torus.Lx - 1), y in 0:(torus.Ly - 1)
        displacement = (x, y)
        point_character = A == identity_matrix ? 1 : parity
        push!(operations, AffineOperation(
            A,
            displacement,
            _affine_permutation(torus, A, displacement),
        ))
        push!(characters,
              _translation_phase(torus, momentum, displacement) * point_character)
    end
    return operations, characters
end

function production_sector_specs(torus::TriangularTorus)
    sectors = NamedTuple[]
    for momentum in production_momentum_representatives(torus)
        star = momentum_star(torus, momentum)
        if length(star) == 1
            for parity in (1, -1)
                push!(sectors, (; momentum, star, parity))
            end
        else
            push!(sectors, (; momentum, star, parity=nothing))
        end
    end
    return sectors
end


function sector_id(momentum::MomentumSector, parity::Union{Nothing,Int}=nothing)
    suffix = parity === nothing ? "" : parity == 1 ? "_c2_p" : "_c2_m"
    return "nx_$(momentum.nx)_ny_$(momentum.ny)$(suffix)"
end
