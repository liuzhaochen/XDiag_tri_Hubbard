@testset "portable JLD2 result records" begin
    record = result_record(
        metadata=Dict("name" => "test"),
        cluster=Dict("N" => 12, "nup" => 6, "ndn" => 6),
        parameters=Dict("U_over_t" => 8.0, "t" => 1.0),
        sector=Dict("nx" => 0, "ny" => 0, "C2_parity" => 1),
        solver=Dict("algorithm" => "test", "residual_norms" => [1e-12, 2e-12]),
        spectrum=Dict(
            "energy" => [-1.0, -0.5],
            "energy_minus_block_ground" => [0.0, 0.5],
            "level_index" => [1, 2],
        ),
    )
    @test validate_result(record)
    mktempdir() do directory
        path = joinpath(directory, "result.jld2")
        @test save_result(path, record) == path
        @test load_result(path) == record
        @test_throws ArgumentError save_result(path, record)
    end
end
