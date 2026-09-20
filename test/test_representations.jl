using XDiag

@testset "momentum stars and XDiag electron blocks" begin
    torus = TriangularTorus(3, 4)
    representatives = production_momentum_representatives(torus)
    @test [(k.nx, k.ny) for k in representatives] ==
          [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (1, 3)]

    expected_stars = [
        [(0, 0)],
        [(0, 1), (0, 3)],
        [(0, 2)],
        [(1, 0), (2, 0)],
        [(1, 1), (2, 3)],
        [(1, 2), (2, 2)],
        [(1, 3), (2, 1)],
    ]
    @test [[(k.nx, k.ny) for k in momentum_star(torus, representative)]
           for representative in representatives] == expected_stars

    translation_dimensions = Dict{Tuple{Int,Int},Int}()
    for momentum in all_momenta(torus)
        block = xdiag_electron_block(torus, 6, 6, momentum)
        translation_dimensions[(momentum.nx, momentum.ny)] = XDiag.dim(block)
    end
    @test sum(values(translation_dimensions)) == binomial(12, 6)^2

    star_weighted_dimension = sum(
        length(momentum_star(torus, momentum)) *
        translation_dimensions[(momentum.nx, momentum.ny)]
        for momentum in representatives
    )
    @test star_weighted_dimension == binomial(12, 6)^2

    for momentum in (MomentumSector(0, 0), MomentumSector(0, 2))
        translation_dimension = translation_dimensions[(momentum.nx, momentum.ny)]
        parity_dimensions = [
            XDiag.dim(xdiag_electron_block(torus, 6, 6, momentum; parity=parity))
            for parity in (-1, 1)
        ]
        @test sum(parity_dimensions) == translation_dimension
    end
end
