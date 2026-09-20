using LinearAlgebra
using XDiag

@testset "direct Hubbard Hamiltonian" begin
    torus = TriangularTorus(3, 4)
    hamiltonian = build_hubbard_hamiltonian(torus; t=1.0, U=8.0)
    @test length(hamiltonian) == 37
    @test isreal(hamiltonian)

    # A zero-hopping block makes the HubbardU convention independently transparent:
    # with one up and one down electron on two sites, double occupancy has energies
    # U,U and separated particles have energies 0,0.
    onsite = XDiag.OpSum(2.5, XDiag.Op("HubbardU"))
    block = XDiag.Electron(2, 1, 1)
    values = sort(real.(LinearAlgebra.eigvals(
        Hermitian(Matrix(XDiag.matrix(onsite, block)))
    )))
    @test isapprox(values, [0.0, 0.0, 2.5, 2.5]; atol=1e-13)

    hopping = XDiag.OpSum(XDiag.Op("Hop", [1, 2]))
    one_electron = XDiag.Electron(2, 1, 0)
    hopping_values = sort(real.(LinearAlgebra.eigvals(
        Hermitian(Matrix(XDiag.matrix(hopping, one_electron)))
    )))
    @test isapprox(hopping_values, [-1.0, 1.0]; atol=1e-13)
end

@testset "full and translation-resolved Hubbard spectra agree" begin
    # This deliberately small torus makes a complete dense reconstruction cheap.
    # The same geometry, hopping convention, fermionic blocks, and translation
    # characters are used by the 3x4 production calculation.
    torus = TriangularTorus(3, 2)
    hamiltonian = build_hubbard_hamiltonian(torus; t=1.0, U=3.7)
    full_block = XDiag.Electron(length(torus), 1, 1)
    full_values = sort(real.(LinearAlgebra.eigvals(Hermitian(Matrix(
        XDiag.matrix(hamiltonian, full_block)
    )))))

    resolved_values = Float64[]
    resolved_dimension = 0
    for momentum in all_momenta(torus)
        block = xdiag_electron_block(torus, 1, 1, momentum)
        resolved_dimension += XDiag.dim(block)
        append!(resolved_values, real.(LinearAlgebra.eigvals(Hermitian(Matrix(
            XDiag.matrix(hamiltonian, block)
        )))))
    end
    sort!(resolved_values)
    @test resolved_dimension == length(full_block) == 36
    @test isapprox(resolved_values, full_values; atol=1e-11)
end
