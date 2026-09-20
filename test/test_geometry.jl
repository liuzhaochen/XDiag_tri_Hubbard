@testset "3x4 triangular torus geometry" begin
    torus = TriangularTorus(3, 4)
    @test length(torus) == 12
    @test [site_index(torus, x, y) for x in 0:2, y in 0:3] ==
          reshape(collect(1:12), 4, 3)'
    @test all(site_coordinate(torus, site_index(torus, x, y)) == (x, y)
              for x in 0:2, y in 0:3)

    bonds = nearest_neighbor_bonds(torus)
    @test length(bonds) == 36
    @test length(Set(bonds)) == 36
    @test all(i < j for (i, j) in bonds)
    @test coordination_numbers(torus, bonds) == fill(6, 12)

    tx = translation_permutation(torus, 1, 0)
    ty = translation_permutation(torus, 0, 1)
    @test tx == [site_index(torus, x + 1, y)
                 for x in 0:2 for y in 0:3]
    @test ty == [site_index(torus, x, y + 1)
                 for x in 0:2 for y in 0:3]
end
