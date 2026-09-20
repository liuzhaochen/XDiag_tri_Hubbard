import XDiag

function xdiag_representation(
    operations::Vector{AffineOperation},
    characters::Vector{ComplexF64},
)
    length(operations) == length(characters) ||
        throw(DimensionMismatch("one character is required per group operation"))
    permutations = [XDiag.Permutation(operation.permutation) for operation in operations]
    group = XDiag.PermutationGroup(permutations)
    return XDiag.Representation(group, characters)
end

function xdiag_electron_block(
    torus::TriangularTorus,
    nup::Int,
    ndn::Int,
    momentum::MomentumSector;
    parity::Union{Nothing,Int}=nothing,
)
    operations, characters = if parity === nothing
        translation_operations_and_characters(torus, momentum)
    else
        little_group_operations_and_characters(torus, momentum, parity)
    end
    representation = xdiag_representation(operations, characters)
    return XDiag.Electron(length(torus), nup, ndn, representation)
end
