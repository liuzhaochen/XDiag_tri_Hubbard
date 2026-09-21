using Dates
using LinearAlgebra
using TOML
using XDiag
using XDiagTriHubbard

function solve_states(operator, block, solver_config)
    requested = Int(solver_config["neigenvalues"])
    residual_limit = Float64(solver_config["residual_acceptance_limit"])
    orthogonality_limit = Float64(solver_config["orthogonality_tolerance"])
    base_seed = Int(solver_config["random_seed"])
    retry_count = Int(solver_config["retry_count"])
    seed_stride = Int(solver_config["retry_seed_stride"])

    solve_seconds = 0.0
    failed_attempts = String[]
    for attempt in 1:retry_count
        seed = base_seed + (attempt - 1) * seed_stride
        println("LOBPCG attempt ", attempt, "/", retry_count, ", seed=", seed)
        flush(stdout)
        timing = @timed XDiag.eigs_lobpcg(
            operator,
            block;
            neigs=requested,
            guard=Int(solver_config["guard_eigenvalues"]),
            tol=Float64(solver_config["precision"]),
            max_iterations=Int(solver_config["max_iterations"]),
            random_seed=seed,
        )
        solve_seconds += timing.time
        candidate = timing.value
        if length(candidate.eigenvalues) < requested ||
           length(candidate.residual_norms) < requested
            push!(failed_attempts, "seed=$seed returned too few states")
            continue
        end
        vectors = XDiag.matrix(candidate.eigenvectors; copy=true)[:, 1:requested]
        values = Float64.(candidate.eigenvalues[1:requested])
        explicit = [
            norm(XDiag.apply(operator, vectors[:, level]) - values[level] * vectors[:, level])
            for level in 1:requested
        ]
        orthogonality = maximum(abs, vectors' * vectors - I)
        accepted = candidate.criterion == "converged" &&
            maximum(candidate.residual_norms[1:requested]) <= residual_limit &&
            maximum(explicit) <= residual_limit &&
            orthogonality <= orthogonality_limit
        if accepted
            return (
                solution=candidate,
                vectors=vectors,
                energies=values,
                explicit_residuals=explicit,
                orthogonality_error=orthogonality,
                seed=seed,
                failed_attempts=failed_attempts,
                solve_seconds=solve_seconds,
            )
        end
        push!(failed_attempts,
              "seed=$seed criterion=$(candidate.criterion) " *
              "reported=$(maximum(candidate.residual_norms[1:requested])) " *
              "explicit=$(maximum(explicit)) orthogonality=$orthogonality")
    end
    error("LOBPCG retries exhausted: $(join(failed_attempts, "; "))")
end

function observable_group(torus, momentum, parity)
    operations, _ = parity === nothing ?
        translation_operations_and_characters(torus, momentum) :
        little_group_operations_and_characters(torus, momentum, parity)
    return xdiag_permutation_group(operations)
end

function stack_level_observables(measurements, key)
    first = measurements[1][key]
    if first isa Number
        return [measurement[key] for measurement in measurements]
    end
    result = Array{eltype(first)}(undef, length(measurements), size(first)...)
    for level in eachindex(measurements)
        selectdim(result, 1, level) .= measurements[level][key]
    end
    return result
end

function git_revision()
    try
        return readchomp(pipeline(`git rev-parse HEAD`; stderr=devnull))
    catch
        return "unavailable"
    end
end

function main(args)
    length(args) == 3 || error(
        "usage: run_observables.jl CONFIG SECTOR_INDEX PARAMETER_INDEX"
    )
    config_path = abspath(args[1])
    sector_index = parse(Int, args[2])
    parameter_index = parse(Int, args[3])
    config = TOML.parsefile(config_path)

    cluster_config = config["cluster"]
    torus = TriangularTorus(Int(cluster_config["Lx"]), Int(cluster_config["Ly"]))
    nup, ndn = Int(cluster_config["nup"]), Int(cluster_config["ndn"])
    nup == ndn == 6 || error("observable production requires Nup=Ndn=6")
    sectors = production_sector_specs(torus)
    length(sectors) == Int(config["sectors"]["expected_count"]) ||
        error("sector count mismatch")
    1 <= sector_index <= length(sectors) || error("sector index out of range")
    sector = sectors[sector_index]
    momentum, parity = sector.momentum, sector.parity
    identifier = sector_id(momentum, parity)

    u_values = Float64.(config["trajectory"]["values"])
    1 <= parameter_index <= length(u_values) || error("parameter index out of range")
    u_value = u_values[parameter_index]
    t_value = Float64(config["trajectory"]["t"])
    u_milli = round(Int, 1000 * u_value)
    filename = "u_$(lpad(u_milli, 5, '0')).jld2"
    output_path = joinpath(abspath(config["output"]["root"]), identifier, filename)
    overwrite = Bool(config["output"]["overwrite"])
    if isfile(output_path) && !overwrite
        println("skip existing: ", output_path)
        return
    end

    println("observable sector=", identifier, ", U/t=", u_value)
    block_timing = @timed xdiag_electron_block(torus, nup, ndn, momentum; parity=parity)
    block = block_timing.value
    hamiltonian = build_hubbard_hamiltonian(torus; t=t_value, U=u_value * t_value)
    operator_timing = @timed XDiag.csr_matrix_32(hamiltonian, block)
    operator = operator_timing.value
    solved = solve_states(operator, block, config["solver"])
    println("energies=", solved.energies)
    flush(stdout)

    reference_path = joinpath(
        abspath(config["observables"]["reference_root"]), identifier, filename
    )
    reference = load_result(reference_path)
    reference_energies = Float64.(reference["spectrum"]["energy"])
    reference_difference = maximum(abs, solved.energies - reference_energies)
    reference_difference <= 1e-7 || error(
        "observable rerun disagrees with coarse spectrum by $reference_difference"
    )

    group = observable_group(torus, momentum, parity)
    measurements = Dict{String,Any}[]
    measurement_seconds = 0.0
    for level in axes(solved.vectors, 2)
        println("measuring level ", level)
        flush(stdout)
        state = XDiag.State(block, solved.vectors[:, level])
        timing = @timed measure_observables(state, torus, group)
        measurement_seconds += timing.time
        push!(measurements, timing.value)
    end

    observable_data = Dict{String,Any}(
        "double_occupancy_per_site" =>
            stack_level_observables(measurements, "double_occupancy_per_site"),
        "spin_displacement_correlation" =>
            stack_level_observables(measurements, "spin_displacement_correlation"),
        "charge_connected_displacement_correlation" => stack_level_observables(
            measurements, "charge_connected_displacement_correlation"
        ),
        "chirality_mean_by_sublattice" =>
            stack_level_observables(measurements, "chirality_mean_by_sublattice"),
        "triangle_sublattice_centers" => measurements[1]["triangle_sublattice_centers"],
        "displacement_convention" => measurements[1]["displacement_convention"],
        "chirality_convention" => measurements[1]["chirality_convention"],
        "measurement_seconds" => measurement_seconds,
        "reference_energy_max_difference" => reference_difference,
    )

    solution = solved.solution
    requested = Int(config["solver"]["neigenvalues"])
    record = result_record(
        metadata=Dict(
            "name" => config["name"],
            "created_utc" => string(now(UTC)),
            "config_path" => relpath(config_path, pwd()),
        ),
        cluster=Dict(
            "Lx" => torus.Lx, "Ly" => torus.Ly, "N" => length(torus),
            "nup" => nup, "ndn" => ndn,
            "bond_count" => length(nearest_neighbor_bonds(torus)),
        ),
        parameters=Dict("U_over_t" => u_value, "t" => t_value, "U" => u_value * t_value),
        sector=Dict(
            "id" => identifier, "nx" => momentum.nx, "ny" => momentum.ny,
            "Ky" => 2pi * momentum.ny / torus.Ly,
            "C2_parity" => parity === nothing ? 0 : parity,
            "star_nx" => [k.nx for k in sector.star],
            "star_ny" => [k.ny for k in sector.star],
            "star_size" => length(sector.star),
            "block_dimension" => XDiag.dim(block),
        ),
        solver=Dict(
            "algorithm" => "eigs_lobpcg", "operator_mode" => "csr32",
            "criterion" => solution.criterion, "iterations" => solution.niterations,
            "requested_eigenvalues" => requested,
            "guard_eigenvalues" => Int(config["solver"]["guard_eigenvalues"]),
            "precision" => Float64(config["solver"]["precision"]),
            "residual_acceptance_limit" => Float64(config["solver"]["residual_acceptance_limit"]),
            "residual_norms" => Float64.(solution.residual_norms[1:requested]),
            "explicit_residual_norms" => solved.explicit_residuals,
            "orthogonality_error" => solved.orthogonality_error,
            "selected_random_seed" => solved.seed,
            "failed_attempts" => solved.failed_attempts,
            "block_construction_seconds" => block_timing.time,
            "operator_construction_seconds" => operator_timing.time,
            "solve_seconds" => solved.solve_seconds,
            "operator_nonzeros" => length(operator.data),
        ),
        spectrum=Dict(
            "energy" => solved.energies,
            "energy_minus_block_ground" => solved.energies .- minimum(solved.energies),
            "level_index" => collect(1:requested),
        ),
        observables=observable_data,
        provenance=Dict(
            "julia_version" => string(VERSION),
            "xdiag_version" => string(pkgversion(XDiag)),
            "git_revision" => git_revision(),
            "slurm_job_id" => get(ENV, "SLURM_JOB_ID", "interactive"),
            "slurm_array_task_id" => get(ENV, "SLURM_ARRAY_TASK_ID", "none"),
        ),
    )
    saved = save_result(output_path, record; overwrite=overwrite)
    println("saved: ", saved)
end

main(ARGS)
