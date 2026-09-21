using JLD2
using TOML
using XDiagTriHubbard

length(ARGS) == 2 || error("usage: summarize_scan.jl CONFIG OUTPUT_JLD2")
config_path = abspath(ARGS[1])
output_path = abspath(ARGS[2])
config = TOML.parsefile(config_path)

cluster = config["cluster"]
torus = TriangularTorus(Int(cluster["Lx"]), Int(cluster["Ly"]))
all_sectors = production_sector_specs(torus)
selected_indices = haskey(config["sectors"], "indices") ?
    Int.(config["sectors"]["indices"]) : collect(eachindex(all_sectors))
sectors = all_sectors[selected_indices]
length(sectors) == Int(config["sectors"]["expected_count"]) ||
    error("selected sector count mismatch")
u_values = Float64.(config["trajectory"]["values"])
nlevels = Int(config["solver"]["neigenvalues"])
raw_root = abspath(config["output"]["root"])

energies = Array{Float64}(undef, length(u_values), length(sectors), nlevels)
reported_residuals = similar(energies)
explicit_residuals = similar(energies)
dimensions = Vector{Int}(undef, length(sectors))
orthogonality = Matrix{Float64}(undef, length(u_values), length(sectors))
iterations = Matrix{Int}(undef, length(u_values), length(sectors))
solve_seconds = Matrix{Float64}(undef, length(u_values), length(sectors))
expected_paths = Set{String}()

for (sector_index, sector) in enumerate(sectors)
    identifier = sector_id(sector.momentum, sector.parity)
    for (u_index, u_value) in enumerate(u_values)
        u_milli = round(Int, 1000 * u_value)
        path = joinpath(raw_root, identifier, "u_$(lpad(u_milli, 5, '0')).jld2")
        push!(expected_paths, abspath(path))
        isfile(path) || error("missing result: $path")
        record = load_result(path)
        record["metadata"]["name"] == config["name"] || error("scan-name mismatch")
        record["cluster"]["N"] == length(torus) || error("cluster mismatch")
        record["sector"]["id"] == identifier || error("sector mismatch")
        isapprox(record["parameters"]["U_over_t"], u_value; atol=1e-12) ||
            error("parameter mismatch")
        record["solver"]["criterion"] == "converged" || error("unconverged result")

        values = Float64.(record["spectrum"]["energy"])
        residual = Float64.(record["solver"]["residual_norms"])
        explicit = Float64.(record["solver"]["explicit_residual_norms"])
        length(values) == length(residual) == length(explicit) == nlevels ||
            error("level-count mismatch")
        energies[u_index, sector_index, :] .= values
        reported_residuals[u_index, sector_index, :] .= residual
        explicit_residuals[u_index, sector_index, :] .= explicit
        dimensions[sector_index] = Int(record["sector"]["block_dimension"])
        orthogonality[u_index, sector_index] = record["solver"]["orthogonality_error"]
        iterations[u_index, sector_index] = Int(record["solver"]["iterations"])
        solve_seconds[u_index, sector_index] = record["solver"]["solve_seconds"]
    end
end

actual_paths = Set(abspath(path) for (root, _, files) in walkdir(raw_root)
                   for file in files if endswith(lowercase(file), ".jld2")
                   for path in (joinpath(root, file),))
unexpected = setdiff(actual_paths, expected_paths)
isempty(unexpected) || error("unexpected result files: $(join(sort!(collect(unexpected)), ", "))")

ground_energy = [minimum(energies[u_index, :, :]) for u_index in eachindex(u_values)]
energy_minus_ground = similar(energies)
for u_index in eachindex(u_values)
    energy_minus_ground[u_index, :, :] .=
        energies[u_index, :, :] .- ground_energy[u_index]
end

summary = Dict{String,Any}(
    "schema_version" => 1,
    "metadata" => Dict(
        "name" => config["name"],
        "config_path" => relpath(config_path, pwd()),
        "record_count" => length(expected_paths),
    ),
    "cluster" => Dict(
        "Lx" => torus.Lx,
        "Ly" => torus.Ly,
        "N" => length(torus),
        "nup" => Int(cluster["nup"]),
        "ndn" => Int(cluster["ndn"]),
    ),
    "grid" => Dict("U_over_t" => u_values),
    "sectors" => Dict(
        "id" => [sector_id(s.momentum, s.parity) for s in sectors],
        "nx" => [s.momentum.nx for s in sectors],
        "ny" => [s.momentum.ny for s in sectors],
        "Ky" => [2pi * s.momentum.ny / torus.Ly for s in sectors],
        "C2_parity" => [s.parity === nothing ? 0 : s.parity for s in sectors],
        "star_size" => [length(s.star) for s in sectors],
        "block_dimension" => dimensions,
    ),
    "spectrum" => Dict(
        "energy" => energies,
        "energy_minus_ground" => energy_minus_ground,
        "ground_energy" => ground_energy,
        "level_index" => collect(1:nlevels),
    ),
    "solver" => Dict(
        "reported_residual_norm" => reported_residuals,
        "explicit_residual_norm" => explicit_residuals,
        "orthogonality_error" => orthogonality,
        "iterations" => iterations,
        "solve_seconds" => solve_seconds,
    ),
)

mkpath(dirname(output_path))
temporary_path = tempname(dirname(output_path)) * ".jld2"
try
    JLD2.jldopen(temporary_path, "w") do file
        file["summary"] = summary
    end
    mv(temporary_path, output_path; force=true)
catch
    isfile(temporary_path) && rm(temporary_path)
    rethrow()
end

println("validated records: ", length(expected_paths))
println("maximum reported residual: ", maximum(reported_residuals))
println("maximum explicit residual: ", maximum(explicit_residuals))
println("maximum orthogonality error: ", maximum(orthogonality))
println("saved summary: ", output_path)

csv_path = replace(output_path, r"\.jld2$" => ".csv")
open(csv_path, "w") do io
    println(io,
        "U_over_t,energy,energy_minus_ground,nx,ny,Ky,C2_parity," *
        "level_index,block_dimension,residual_norm,explicit_residual_norm," *
        "orthogonality_error"
    )
    for u_index in eachindex(u_values), sector_index in eachindex(sectors), level in 1:nlevels
        sector = sectors[sector_index]
        parity = sector.parity === nothing ? 0 : sector.parity
        println(io, join((
            u_values[u_index],
            energies[u_index, sector_index, level],
            energy_minus_ground[u_index, sector_index, level],
            sector.momentum.nx,
            sector.momentum.ny,
            2pi * sector.momentum.ny / torus.Ly,
            parity,
            level,
            dimensions[sector_index],
            reported_residuals[u_index, sector_index, level],
            explicit_residuals[u_index, sector_index, level],
            orthogonality[u_index, sector_index],
        ), ','))
    end
end
println("saved long-form table: ", csv_path)
