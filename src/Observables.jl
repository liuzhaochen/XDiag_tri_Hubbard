import LinearAlgebra
import XDiag

"""Counterclockwise elementary triangles and centers in `(a1,a2)` coordinates."""
function oriented_triangles(torus::TriangularTorus)
    triangles = NamedTuple[]
    for x in 0:(torus.Lx - 1), y in 0:(torus.Ly - 1)
        push!(triangles, (
            sites=(site_index(torus, x, y), site_index(torus, x + 1, y),
                   site_index(torus, x, y + 1)),
            center=(x + 1 / 3, y + 1 / 3),
            orientation="up",
            cell=(x, y),
        ))
        push!(triangles, (
            sites=(site_index(torus, x + 1, y), site_index(torus, x + 1, y + 1),
                   site_index(torus, x, y + 1)),
            center=(x + 2 / 3, y + 2 / 3),
            orientation="down",
            cell=(x, y),
        ))
    end
    return triangles
end

site_fractional_coordinates(torus::TriangularTorus) =
    [Float64.(collect(site_coordinate(torus, site))) for site in 1:length(torus)]

"""Fourier transform a finite correlation matrix, normalized per object."""
function structure_factor(correlation::AbstractMatrix, coordinates, q::NTuple{2,<:Real})
    nobjects = length(coordinates)
    size(correlation) == (nobjects, nobjects) ||
        throw(DimensionMismatch("correlation matrix and coordinates disagree"))
    phases = ComplexF64[
        cis(2pi * (q[1] * coordinate[1] + q[2] * coordinate[2]))
        for coordinate in coordinates
    ]
    return LinearAlgebra.dot(phases, correlation * phases) / nobjects
end

"""Fourier transform a translation-averaged site correlation `C(dx,dy)`."""
function displacement_structure_factor(correlation::AbstractMatrix, q::NTuple{2,<:Real})
    Lx, Ly = size(correlation)
    return sum(
        cis(2pi * (q[1] * dx + q[2] * dy)) * correlation[dx + 1, dy + 1]
        for dx in 0:(Lx - 1), dy in 0:(Ly - 1)
    )
end

"""Fourier transform two-sublattice triangle correlations, per triangle."""
function triangle_structure_factor(
    correlation::AbstractArray{<:Number,4}, centers, q::NTuple{2,<:Real},
)
    size(correlation, 1) == size(correlation, 2) == length(centers) ||
        throw(DimensionMismatch("triangle sublattices disagree"))
    Lx, Ly = size(correlation, 3), size(correlation, 4)
    value = zero(ComplexF64)
    for first in eachindex(centers), second in eachindex(centers)
        center_shift = (centers[second][1] - centers[first][1],
                        centers[second][2] - centers[first][2])
        for dx in 0:(Lx - 1), dy in 0:(Ly - 1)
            phase = cis(2pi * (q[1] * (dx + center_shift[1]) +
                               q[2] * (dy + center_shift[2])))
            value += phase * correlation[first, second, dx + 1, dy + 1]
        end
    end
    return value / length(centers)
end

function _density_operator(site::Int)
    return XDiag.OpSum(XDiag.Op("Nup", site)) + XDiag.Op("Ndn", site)
end

function _spin_square_operator(site::Int)
    return 0.75 * (XDiag.OpSum(XDiag.Op("Nup", site)) + XDiag.Op("Ndn", site)) -
           1.5 * XDiag.Op("Nupdn", site)
end

_as_opsum(operator::XDiag.OpSum) = operator
_as_opsum(operator) = XDiag.OpSum(operator)
_invariant(operator, group) = XDiag.symmetrize(_as_opsum(operator), group)

function _expect_invariant(state, operator, group)
    return XDiag.inner(_invariant(operator, group), state)
end

function _site_correlations(state, torus::TriangularTorus, group)
    spin = Matrix{ComplexF64}(undef, torus.Lx, torus.Ly)
    charge = similar(spin)
    origin = site_index(torus, 0, 0)
    density_origin = _density_operator(origin)
    for dx in 0:(torus.Lx - 1), dy in 0:(torus.Ly - 1)
        destination = site_index(torus, dx, dy)
        spin_operator = destination == origin ?
            _spin_square_operator(origin) : XDiag.Op("SdotS", [origin, destination])
        density_destination = _density_operator(destination)
        spin[dx + 1, dy + 1] = _expect_invariant(state, spin_operator, group)
        charge[dx + 1, dy + 1] =
            _expect_invariant(state, density_origin * density_destination, group) - 1
    end
    return spin, charge
end

function _reference_triangles(torus::TriangularTorus)
    triangles = oriented_triangles(torus)
    return [only(triangle for triangle in triangles
                 if triangle.cell == (0, 0) && triangle.orientation == orientation)
            for orientation in ("up", "down")]
end

function _chirality_operator(sites)
    i, j, k = sites
    first_bond = XDiag.Op("SdotS", [i, j])
    second_bond = XDiag.Op("SdotS", [j, k])
    # [Si.Sj, Sj.Sk] = -i Si.(Sj x Sk).
    return 1.0im * (first_bond * second_bond - second_bond * first_bond)
end

function _chirality_means(state, torus::TriangularTorus, group)
    references = _reference_triangles(torus)
    centers = [collect(reference.center) for reference in references]
    means = ComplexF64[
        _expect_invariant(state, _chirality_operator(reference.sites), group)
        for reference in references
    ]
    return means, centers
end

"""Measure translation-averaged observables in a symmetry-adapted state.

Site correlations are indexed by `(dx,dy)`. The two uniform chirality means are
for the up/down elementary triangles, both oriented counterclockwise in Cartesian
space.
"""
function measure_observables(state, torus::TriangularTorus, group)
    double_occupancy = real(XDiag.inner(XDiag.Op("HubbardU"), state)) / length(torus)
    spin, charge = _site_correlations(state, torus, group)
    chirality_mean, triangle_centers = _chirality_means(state, torus, group)
    return Dict{String,Any}(
        "double_occupancy_per_site" => double_occupancy,
        "spin_displacement_correlation" => spin,
        "charge_connected_displacement_correlation" => charge,
        "chirality_mean_by_sublattice" => chirality_mean,
        "triangle_sublattice_centers" => triangle_centers,
        "displacement_convention" => "dx=0:Lx-1, dy=0:Ly-1 in (a1,a2)",
        "chirality_convention" =>
            "S_i dot (S_j cross S_k), sites counterclockwise in Cartesian space",
    )
end

function chirality_fourier_operator(
    torus::TriangularTorus, q::NTuple{2,<:Real},
)
    triangles = oriented_triangles(torus)
    operator = XDiag.OpSum()
    normalization = sqrt(length(triangles))
    for triangle in triangles
        phase = cis(2pi * (q[1] * triangle.center[1] + q[2] * triangle.center[2]))
        operator += (phase / normalization) * _chirality_operator(triangle.sites)
    end
    return operator
end

chirality_fourier_operator(torus::TriangularTorus, momentum::MomentumSector) =
    chirality_fourier_operator(
        torus, (momentum.nx / torus.Lx, momentum.ny / torus.Ly)
    )

function _chirality_sublattice_fourier_operator(
    torus::TriangularTorus, momentum::MomentumSector, orientation::String,
)
    q = (momentum.nx / torus.Lx, momentum.ny / torus.Ly)
    triangles = filter(t -> t.orientation == orientation, oriented_triangles(torus))
    operator = XDiag.OpSum()
    normalization = sqrt(length(torus))
    for triangle in triangles
        x, y = triangle.cell
        phase = cis(2pi * (q[1] * x + q[2] * y))
        operator += (phase / normalization) * _chirality_operator(triangle.sites)
    end
    return operator
end

"""Scalar-chirality structure factor at all allowed cluster momenta.

The input block must use translation symmetry only. `chi_q` is normalized by
`1/sqrt(Ntriangles)`, so `structure_factor = <chi_q^dagger chi_q>` is per triangle.
"""
function measure_chirality_structure(
    vector::AbstractVector,
    input_block,
    torus::TriangularTorus,
    input_momentum::MomentumSector,
)
    normalized = ComplexF64.(vector) ./ LinearAlgebra.norm(vector)
    full = Matrix{Float64}(undef, torus.Lx, torus.Ly)
    connected = similar(full)
    sublattice_gram = Array{ComplexF64}(undef, 2, 2, torus.Lx, torus.Ly)
    mean_per_triangle = 0.0 + 0.0im
    ntriangles = 2 * length(torus)
    orientations = ("up", "down")
    offsets = ((1 / 3, 1 / 3), (2 / 3, 2 / 3))
    for q_momentum in all_momenta(torus)
        output_momentum = MomentumSector(
            mod(input_momentum.nx - q_momentum.nx, torus.Lx),
            mod(input_momentum.ny - q_momentum.ny, torus.Ly),
        )
        output_block = xdiag_electron_block(torus, 6, 6, output_momentum)
        outputs = Vector{ComplexF64}[]
        for orientation in orientations
            operator = _chirality_sublattice_fourier_operator(
                torus, q_momentum, orientation
            )
            output = zeros(ComplexF64, XDiag.dim(output_block))
            XDiag.apply(operator, input_block, normalized, output_block, output)
            push!(outputs, output)
        end
        gram = ComplexF64[
            LinearAlgebra.dot(outputs[first], outputs[second])
            for first in 1:2, second in 1:2
        ]
        sublattice_gram[:, :, q_momentum.nx + 1, q_momentum.ny + 1] .= gram
        q = (q_momentum.nx / torus.Lx, q_momentum.ny / torus.Ly)
        weights = ComplexF64[
            cis(2pi * (q[1] * offset[1] + q[2] * offset[2])) / sqrt(2)
            for offset in offsets
        ]
        total_output = weights[1] * outputs[1] + weights[2] * outputs[2]
        value = real(LinearAlgebra.dot(total_output, total_output))
        full[q_momentum.nx + 1, q_momentum.ny + 1] = value
        if q_momentum == MomentumSector(0, 0)
            mean_per_triangle =
                LinearAlgebra.dot(normalized, total_output) / sqrt(ntriangles)
            connected[1, 1] = value - ntriangles * abs2(mean_per_triangle)
        else
            connected[q_momentum.nx + 1, q_momentum.ny + 1] = value
        end
    end
    return Dict{String,Any}(
        "chirality_mean_per_triangle" => mean_per_triangle,
        "chirality_structure_factor" => full,
        "chirality_connected_structure_factor" => connected,
        "chirality_sublattice_gram" => sublattice_gram,
        "triangle_sublattice_offsets" => [collect(offset) for offset in offsets],
        "momentum_nx" => collect(0:(torus.Lx - 1)),
        "momentum_ny" => collect(0:(torus.Ly - 1)),
        "normalization" =>
            "chi_(a,q)=sum_cell exp(i q.R) chi_(a,R)/sqrt(N); " *
            "chi_q=sum_a exp(i q.d_a) chi_(a,q)/sqrt(2)",
    )
end
