using JLD2
using LinearAlgebra
using TOML
using XDiagTriHubbard

function real_checked(value, label; tolerance=1e-7)
    abs(imag(value)) <= tolerance || error("complex $label: $value")
    return real(value)
end

function high_symmetry_momenta(torus)
    return [
        ("Gamma", (0.0, 0.0)),
        ("Kcl", (-1 / 3, 1 / 4)),      # allowed (2,1), closest to K=(-1/3,1/3)
        ("M1cl", (1 / 3, 0.0)),        # allowed (1,0), closest to M1=(1/2,0)
        ("M2", (0.0, -1 / 2)),         # allowed (0,2), exact M2
        ("M3cl", (1 / 3, 1 / 2)),      # allowed (1,2), closest to M3=(1/2,1/2)
    ]
end

reciprocal_norm2(q1, q2) = q1^2 - q1 * q2 + q2^2

function folded_fractional_momentum(torus, nx, ny)
    candidates = [
        (nx / torus.Lx + shift_x, ny / torus.Ly + shift_y)
        for shift_x in -1:1, shift_y in -1:1
    ]
    return candidates[argmin(reciprocal_norm2(q...) for q in candidates)]
end

function chirality_from_sublattice_gram(gram, torus, q)
    nx = mod(round(Int, q[1] * torus.Lx), torus.Lx)
    ny = mod(round(Int, q[2] * torus.Ly), torus.Ly)
    offsets = ((1 / 3, 1 / 3), (2 / 3, 2 / 3))
    weights = ComplexF64[
        cis(2pi * (q[1] * offset[1] + q[2] * offset[2])) / sqrt(2)
        for offset in offsets
    ]
    matrix = gram[:, :, nx + 1, ny + 1]
    return real_checked(LinearAlgebra.dot(weights, matrix * weights), "chirality SF")
end

function main(args)
    length(args) == 2 || error("usage: summarize_observables.jl CONFIG OUTPUT_JLD2")
    config_path, output_path = abspath(args[1]), abspath(args[2])
    config = TOML.parsefile(config_path)
    torus = TriangularTorus(3, 4)
    sectors = production_sector_specs(torus)
    u_values = Float64.(config["trajectory"]["values"])
    nlevels = Int(config["solver"]["neigenvalues"])
    raw_root = abspath(config["output"]["root"])
    high_q = high_symmetry_momenta(torus)

    records = Dict{Tuple{Int,Int},Any}()
    expected_paths = Set{String}()
    max_explicit_residual = 0.0
    max_orthogonality_error = 0.0
    max_reference_difference = 0.0
    retried_records = 0
    for (sector_index, sector) in enumerate(sectors), (u_index, u_value) in enumerate(u_values)
        identifier = sector_id(sector.momentum, sector.parity)
        path = joinpath(raw_root, identifier,
                        "u_$(lpad(round(Int, 1000 * u_value), 5, '0')).jld2")
        push!(expected_paths, abspath(path))
        isfile(path) || error("missing observable result: $path")
        record = load_result(path)
        record["metadata"]["name"] == config["name"] || error("scan-name mismatch")
        record["solver"]["criterion"] == "converged" || error("unconverged record")
        residual = maximum(record["solver"]["explicit_residual_norms"])
        residual <= 1e-8 ||
            error("residual gate failed: $path")
        reference_difference = record["observables"]["reference_energy_max_difference"]
        reference_difference <= 1e-7 ||
            error("reference spectrum mismatch: $path")
        max_explicit_residual = max(max_explicit_residual, residual)
        max_orthogonality_error = max(
            max_orthogonality_error, record["solver"]["orthogonality_error"]
        )
        max_reference_difference = max(max_reference_difference, reference_difference)
        retried_records += !isempty(record["solver"]["failed_attempts"])
        records[(u_index, sector_index)] = record
    end

    ground_sector_indices = Int.(config["observables"]["ground_sector_indices"])
    ground_energies = Float64[
        minimum(
            minimum(Float64.(records[(u_index, sector_index)]["spectrum"]["energy"]))
            for sector_index in eachindex(sectors)
        )
        for u_index in eachindex(u_values)
    ]

    state_csv = replace(output_path, r"\.jld2$" => "_states.csv")
    minimum_allowed_spin = Inf
    minimum_allowed_charge = Inf
    maximum_abs_chirality_mean = 0.0
    mkpath(dirname(output_path))
    open(state_csv, "w") do io
        columns = [
            "U_over_t", "energy", "energy_minus_ground", "nx", "ny", "C2_parity",
            "sector_id", "level_index", "double_occupancy", "chirality_mean",
        ]
        for (name, _) in high_q
            append!(columns, ["spin_$name", "charge_$name"])
        end
        println(io, join(columns, ','))
        for (u_index, u_value) in enumerate(u_values), sector_index in eachindex(sectors)
            record = records[(u_index, sector_index)]
            sector = sectors[sector_index]
            observables = record["observables"]
            energies = Float64.(record["spectrum"]["energy"])
            for level in 1:nlevels
                spin = observables["spin_displacement_correlation"][level, :, :]
                charge = observables["charge_connected_displacement_correlation"][level, :, :]
                chi_sublattice = observables["chirality_mean_by_sublattice"][level, :]
                chi_mean = real_checked(sum(chi_sublattice) / 2, "chirality mean")
                maximum_abs_chirality_mean = max(maximum_abs_chirality_mean, abs(chi_mean))
                for nx in 0:2, ny in 0:3
                    q = (nx / 3, ny / 4)
                    spin_allowed = real_checked(
                        displacement_structure_factor(spin, q), "allowed spin"
                    )
                    charge_allowed = real_checked(
                        displacement_structure_factor(charge, q), "allowed charge"
                    )
                    minimum_allowed_spin = min(minimum_allowed_spin, spin_allowed)
                    minimum_allowed_charge = min(minimum_allowed_charge, charge_allowed)
                end
                row = Any[
                    u_value, energies[level], energies[level] - ground_energies[u_index],
                    sector.momentum.nx, sector.momentum.ny,
                    sector.parity === nothing ? 0 : sector.parity,
                    sector_id(sector.momentum, sector.parity), level,
                    observables["double_occupancy_per_site"][level], chi_mean,
                ]
                for (name, q) in high_q
                    spin_value = real_checked(
                        displacement_structure_factor(spin, q), "spin $name"
                    )
                    charge_value = real_checked(
                        displacement_structure_factor(charge, q), "charge $name"
                    )
                    spin_value >= -1e-7 || error("negative spin structure factor")
                    charge_value >= -1e-7 || error("negative charge structure factor")
                    append!(row, [spin_value, charge_value])
                end
                println(io, join(row, ','))
            end
        end
    end

    ground_chirality = Dict{Int,Any}()
    ground_csv = replace(output_path, r"\.jld2$" => "_ground.csv")
    chirality_csv = replace(output_path, r"\.jld2$" => "_chirality_grid.csv")
    chirality_gram_csv = replace(
        output_path, r"\.jld2$" => "_chirality_sublattice_gram.csv"
    )
    open(ground_csv, "w") do ground_io
        header = ["U_over_t", "energy", "nx", "ny", "C2_parity", "double_occupancy",
                  "chirality_mean", "chirality_projected_abs", "chirality_sf_Gamma",
                  "chirality_sf_Kcl", "chirality_sf_M1cl", "chirality_sf_M2",
                  "chirality_sf_M3cl", "chirality_sf_max"]
        for (name, _) in high_q
            append!(header, ["spin_$name", "charge_$name"])
        end
        println(ground_io, join(header, ','))
        open(chirality_csv, "w") do chi_io
            println(chi_io, "U_over_t,ground_nx,ground_ny,q_nx,q_ny,q1,q2,full,connected")
            for (u_index, u_value) in enumerate(u_values)
                sector_index = ground_sector_indices[u_index]
                sector = sectors[sector_index]
                record = records[(u_index, sector_index)]
                observables = record["observables"]
                spin = observables["spin_displacement_correlation"][1, :, :]
                charge = observables["charge_connected_displacement_correlation"][1, :, :]
                chi_path = joinpath(raw_root, "ground_chirality",
                    "u_$(lpad(round(Int, 1000 * u_value), 5, '0')).jld2")
                push!(expected_paths, abspath(chi_path))
                isfile(chi_path) || error("missing ground chirality result: $chi_path")
                chi_record = load_result(chi_path)
                chi_residual = chi_record["solver"]["explicit_residual_norms"][1]
                chi_residual <= 1e-8 ||
                    error("ground chirality residual failed")
                chi_reference_difference =
                    chi_record["observables"]["reference_energy_difference"]
                chi_reference_difference <= 1e-7 ||
                    error("ground chirality reference mismatch")
                max_explicit_residual = max(max_explicit_residual, chi_residual)
                max_reference_difference = max(
                    max_reference_difference, chi_reference_difference
                )
                chi = chi_record["observables"]["chirality"]
                gram = ComplexF64.(chi["chirality_sublattice_gram"])
                for nx in 0:2, ny in 0:3
                    matrix = gram[:, :, nx + 1, ny + 1]
                    isapprox(matrix, matrix'; atol=1e-8) ||
                        error("non-Hermitian chirality sublattice Gram matrix")
                    minimum(LinearAlgebra.eigvals(LinearAlgebra.Hermitian(matrix))) >= -1e-8 ||
                        error("non-positive chirality sublattice Gram matrix")
                end
                first_bz = Matrix{Float64}(undef, torus.Lx, torus.Ly)
                connected = similar(first_bz)
                for nx in 0:2, ny in 0:3
                    q = folded_fractional_momentum(torus, nx, ny)
                    first_bz[nx + 1, ny + 1] =
                        chirality_from_sublattice_gram(gram, torus, q)
                    connected[nx + 1, ny + 1] = first_bz[nx + 1, ny + 1]
                end
                mean_per_triangle = real_checked(
                    chi["chirality_mean_per_triangle"], "ground chi"
                )
                connected[1, 1] -= 2 * length(torus) * abs2(mean_per_triangle)
                minimum(first_bz) >= -1e-8 || error("negative chirality structure factor")
                minimum(connected) >= -1e-8 || error("negative connected chirality SF")
                ground_chirality[u_index] = chi_record
                for nx in 0:2, ny in 0:3
                    q = folded_fractional_momentum(torus, nx, ny)
                    println(chi_io, join((
                        u_value, sector.momentum.nx, sector.momentum.ny,
                        nx, ny, q[1], q[2],
                        first_bz[nx + 1, ny + 1], connected[nx + 1, ny + 1],
                    ), ','))
                end
                chi_values = [chirality_from_sublattice_gram(gram, torus, q)
                              for (_, q) in high_q]
                chi_values[1] -= 2 * length(torus) * abs2(mean_per_triangle)
                projected = maximum(abs, chi_record["observables"][
                    "projected_uniform_chirality_eigenvalues"
                ])
                row = Any[
                    u_value, ground_energies[u_index], sector.momentum.nx,
                    sector.momentum.ny, sector.parity === nothing ? 0 : sector.parity,
                    observables["double_occupancy_per_site"][1],
                    mean_per_triangle, projected,
                    chi_values..., maximum(connected),
                ]
                for (_, q) in high_q
                    append!(row, [
                        real_checked(displacement_structure_factor(spin, q), "ground spin"),
                        real_checked(displacement_structure_factor(charge, q), "ground charge"),
                    ])
                end
                println(ground_io, join(row, ','))
            end
        end
    end

    open(chirality_gram_csv, "w") do io
        println(io, "tag,U_over_t,sublattice_a,sublattice_b,q_nx,q_ny,real,imag")
        for (tag, u_index) in (("U8_k10", 5), ("U14_Gamma", 17))
            gram = ComplexF64.(ground_chirality[u_index]["observables"]["chirality"][
                "chirality_sublattice_gram"
            ])
            for first in 1:2, second in 1:2, nx in 0:2, ny in 0:3
                value = gram[first, second, nx + 1, ny + 1]
                println(io, join((
                    tag, u_values[u_index], first, second, nx, ny,
                    real(value), imag(value),
                ), ','))
            end
        end
    end

    heatmap_csv = replace(output_path, r"\.jld2$" => "_heatmap_correlations.csv")
    open(heatmap_csv, "w") do io
        println(io, "tag,U_over_t,nx,ny,C2_parity,observable,dx,dy,real,imag")
        selections = [
            ("U8_k10", 5, 6),
            ("U14_Gamma", 17, 1),
        ]
        for (tag, u_index, sector_index) in selections
            record = records[(u_index, sector_index)]
            sector = sectors[sector_index]
            for (observable_name, key) in (
                ("spin", "spin_displacement_correlation"),
                ("charge", "charge_connected_displacement_correlation"),
            )
                correlation = record["observables"][key][1, :, :]
                for dx in 0:2, dy in 0:3
                    value = correlation[dx + 1, dy + 1]
                    println(io, join((
                        tag, u_values[u_index], sector.momentum.nx, sector.momentum.ny,
                        sector.parity === nothing ? 0 : sector.parity,
                        observable_name, dx, dy, real(value), imag(value),
                    ), ','))
                end
            end
        end
    end

    actual_paths = Set(abspath(path) for (root, _, files) in walkdir(raw_root)
                       for file in files if endswith(lowercase(file), ".jld2")
                       for path in (joinpath(root, file),))
    unexpected = setdiff(actual_paths, expected_paths)
    isempty(unexpected) || error("unexpected observable files: $(join(sort!(collect(unexpected)), ", "))")

    summary = Dict{String,Any}(
        "metadata" => Dict(
            "name" => config["name"], "record_count" => length(expected_paths),
            "state_csv" => basename(state_csv), "ground_csv" => basename(ground_csv),
            "chirality_csv" => basename(chirality_csv), "heatmap_csv" => basename(heatmap_csv),
            "chirality_gram_csv" => basename(chirality_gram_csv),
        ),
        "grid" => Dict("U_over_t" => u_values),
        "ground" => Dict(
            "energy" => ground_energies, "sector_index" => ground_sector_indices,
        ),
        "validation" => Dict(
            "max_explicit_residual" => max_explicit_residual,
            "max_orthogonality_error" => max_orthogonality_error,
            "max_reference_energy_difference" => max_reference_difference,
            "retried_sector_records" => retried_records,
            "minimum_allowed_spin_structure_factor" => minimum_allowed_spin,
            "minimum_allowed_charge_structure_factor" => minimum_allowed_charge,
            "maximum_abs_uniform_chirality_mean" => maximum_abs_chirality_mean,
        ),
    )
    JLD2.jldopen(output_path, "w") do file
        file["summary"] = summary
    end
    println("validated observable records: ", length(expected_paths))
    println("max explicit residual: ", max_explicit_residual)
    println("max orthogonality error: ", max_orthogonality_error)
    println("max reference energy difference: ", max_reference_difference)
    println("retried sector records: ", retried_records)
    println("minimum allowed spin SF: ", minimum_allowed_spin)
    println("minimum allowed charge SF: ", minimum_allowed_charge)
    println("max |uniform chirality| across retained states: ", maximum_abs_chirality_mean)
    println("saved: ", output_path)
    println("saved: ", state_csv)
    println("saved: ", ground_csv)
    println("saved: ", chirality_csv)
    println("saved: ", chirality_gram_csv)
    println("saved: ", heatmap_csv)
end

main(ARGS)
