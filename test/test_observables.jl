using LinearAlgebra
using XDiag

@testset "observable geometry and Fourier conventions" begin
    torus = TriangularTorus(3, 4)
    triangles = oriented_triangles(torus)
    @test length(triangles) == 24
    @test count(t -> t.orientation == "up", triangles) == 12
    @test count(t -> t.orientation == "down", triangles) == 12
    @test length(unique(Tuple(sort(collect(t.sites))) for t in triangles)) == 24
    @test all(length(unique(t.sites)) == 3 for t in triangles)

    coordinates = site_fractional_coordinates(torus)
    identity_correlation = Matrix{Float64}(I, 12, 12)
    @test isapprox(structure_factor(identity_correlation, coordinates, (0.0, 0.0)), 1)
    @test isapprox(structure_factor(identity_correlation, coordinates, (2 / 3, 1 / 4)), 1)

    uniform = ones(ComplexF64, 12, 12)
    @test isapprox(structure_factor(uniform, coordinates, (0.0, 0.0)), 12)
    @test abs(structure_factor(uniform, coordinates, (1 / 3, 0.0))) <= 1e-12
end

@testset "scalar-chirality operator conventions" begin
    torus = TriangularTorus(3, 2)
    block = XDiag.Electron(length(torus), 3, 3)
    triangle = first(oriented_triangles(torus))
    local_chi = XDiagTriHubbard._chirality_operator(triangle.sites)
    local_matrix = Matrix(XDiag.matrix(local_chi, block))
    @test isapprox(local_matrix, local_matrix'; atol=1e-13)

    q = (1 / 3, 1 / 2)
    minus_q = (-1 / 3, -1 / 2)
    chi_q = Matrix(XDiag.matrix(chirality_fourier_operator(torus, q), block))
    chi_minus_q = Matrix(XDiag.matrix(chirality_fourier_operator(torus, minus_q), block))
    @test isapprox(chi_q', chi_minus_q; atol=1e-13)

    # A real normalized vector is time-reversal invariant in the occupation basis.
    vector = normalize!(collect(1.0:XDiag.dim(block)))
    @test isapprox(norm(chi_q * vector)^2, norm(chi_minus_q * vector)^2; atol=1e-12)
end
