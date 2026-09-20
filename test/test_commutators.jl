using LinearAlgebra
using XDiag

function fixed_masks(nsites::Int, number::Int)
    return [mask for mask in 0:((Int(1) << nsites) - 1)
            if count_ones(mask) == number]
end

function hopping_sign(mask::Int, destination::Int, origin::Int)
    lower, upper = min(destination, origin), max(destination, origin)
    between = 0
    for site in (lower + 1):(upper - 1)
        between += (mask >> (site - 1)) & 1
    end
    return iseven(between) ? 1.0 : -1.0
end

function permute_mask_with_sign(mask::Int, permutation::Vector{Int})
    occupied = [site for site in eachindex(permutation)
                if ((mask >> (site - 1)) & 1) == 1]
    images = permutation[occupied]
    inversions = 0
    for i in eachindex(images), j in (i + 1):length(images)
        inversions += images[i] > images[j]
    end
    mapped = 0
    for site in images
        mapped |= Int(1) << (site - 1)
    end
    return mapped, iseven(inversions) ? 1.0 : -1.0
end

function independent_hubbard_matrix(
    torus::TriangularTorus,
    nup::Int,
    ndn::Int;
    t::Float64,
    U::Float64,
)
    up_masks, dn_masks = fixed_masks(length(torus), nup), fixed_masks(length(torus), ndn)
    basis = [(up, dn) for up in up_masks for dn in dn_masks]
    index = Dict(state => position for (position, state) in enumerate(basis))
    matrix = zeros(Float64, length(basis), length(basis))
    for (column, (up, dn)) in enumerate(basis)
        matrix[column, column] += U * count_ones(up & dn)
        for (i, j) in nearest_neighbor_bonds(torus), spin in 1:2
            mask = spin == 1 ? up : dn
            for (destination, origin) in ((i, j), (j, i))
                destination_bit = Int(1) << (destination - 1)
                origin_bit = Int(1) << (origin - 1)
                mask & origin_bit == 0 && continue
                mask & destination_bit != 0 && continue
                moved = xor(xor(mask, origin_bit), destination_bit)
                new_state = spin == 1 ? (moved, dn) : (up, moved)
                row = index[new_state]
                matrix[row, column] += -t * hopping_sign(mask, destination, origin)
            end
        end
    end
    return matrix, basis, index
end

function fermionic_permutation_matrix(
    basis::Vector{Tuple{Int,Int}},
    index::Dict{Tuple{Int,Int},Int},
    permutation::Vector{Int},
)
    matrix = zeros(Float64, length(basis), length(basis))
    for (column, (up, dn)) in enumerate(basis)
        mapped_up, sign_up = permute_mask_with_sign(up, permutation)
        mapped_dn, sign_dn = permute_mask_with_sign(dn, permutation)
        matrix[index[(mapped_up, mapped_dn)], column] = sign_up * sign_dn
    end
    return matrix
end

@testset "independent fermionic symmetry commutators" begin
    torus = TriangularTorus(3, 2)
    U = 3.7
    independent, basis, index = independent_hubbard_matrix(
        torus, 1, 1; t=1.0, U=U
    )
    @test LinearAlgebra.ishermitian(independent)

    inversion = only(A for A in surviving_point_group(torus) if A == [-1 0; 0 -1])
    permutations = (
        translation_permutation(torus, 1, 0),
        translation_permutation(torus, 0, 1),
        point_permutation(torus, inversion),
    )
    for permutation in permutations
        symmetry = fermionic_permutation_matrix(basis, index, permutation)
        @test isapprox(symmetry' * symmetry, I; atol=1e-13)
        @test maximum(abs, independent * symmetry - symmetry * independent) <= 1e-13
    end

    xdiag = build_hubbard_hamiltonian(torus; t=1.0, U=U)
    block = XDiag.Electron(length(torus), 1, 1)
    xdiag_values = sort(real.(LinearAlgebra.eigvals(Hermitian(Matrix(
        XDiag.matrix(xdiag, block)
    )))))
    independent_values = sort(LinearAlgebra.eigvals(Hermitian(independent)))
    @test isapprox(xdiag_values, independent_values; atol=1e-11)
end
