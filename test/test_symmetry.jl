@testset "triangular point group and surviving torus symmetry" begin
    torus = TriangularTorus(3, 4)
    d6 = integer_point_group()
    @test length(d6) == 12

    point_group = surviving_point_group(torus)
    @test length(point_group) == 2
    @test Set(Tuple(vec(A)) for A in point_group) ==
          Set((Tuple(vec([1 0; 0 1])), Tuple(vec([-1 0; 0 -1]))))
    @test all(preserves_bonds(torus, A) for A in point_group)

    point_permutations = [point_permutation(torus, A) for A in point_group]
    @test is_closed_permutation_group(point_permutations)

    translations = [translation_permutation(torus, x, y)
                    for x in 0:2 for y in 0:3]
    @test length(Set(Tuple(p) for p in translations)) == 12
    @test is_closed_permutation_group(translations)

    compose(left, right) = [left[right[site]] for site in eachindex(right)]
    identity = collect(1:length(torus))
    tx = translation_permutation(torus, 1, 0)
    ty = translation_permutation(torus, 0, 1)
    c2 = point_permutation(torus, [-1 0; 0 -1])
    @test compose(tx, compose(tx, tx)) == identity
    @test compose(ty, compose(ty, compose(ty, ty))) == identity
    @test compose(c2, c2) == identity
    @test compose(tx, ty) == compose(ty, tx)
    @test compose(c2, compose(tx, c2)) == translation_permutation(torus, -1, 0)
    @test compose(c2, compose(ty, c2)) == translation_permutation(torus, 0, -1)
end
