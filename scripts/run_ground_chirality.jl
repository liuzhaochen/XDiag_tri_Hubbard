using Dates
using LinearAlgebra
using TOML
using XDiag
using XDiagTriHubbard

function trusted_ground(operator, block, solver_config)
    for attempt in 1:Int(solver_config["retry_count"])
        seed = Int(solver_config["random_seed"]) +
               (attempt - 1) * Int(solver_config["retry_seed_stride"])
        timing = @timed XDiag.eigs_lobpcg(
            operator, block; neigs=1, guard=2,
            tol=Float64(solver_config["precision"]),
            max_iterations=Int(solver_config["max_iterations"]), random_seed=seed,
        )
        solution = timing.value
        if solution.criterion != "converged" || isempty(solution.eigenvalues)
            continue
        end
        vector = XDiag.matrix(solution.eigenvectors; copy=true)[:, 1]
        energy = Float64(solution.eigenvalues[1])
        residual = norm(XDiag.apply(operator, vector) - energy * vector)
        if residual <= Float64(solver_config["residual_acceptance_limit"])
            return (; energy, vector, residual, solution, seed, seconds=timing.time)
        end
    end
    error("ground-state LOBPCG retries exhausted")
end

function main(args)
    length(args) == 2 || error(
        "usage: run_ground_chirality.jl CONFIG PARAMETER_INDEX"
    )
    config_path = abspath(args[1])
    parameter_index = parse(Int, args[2])
    config = TOML.parsefile(config_path)
    torus = TriangularTorus(3, 4)
    u_values = Float64.(config["trajectory"]["values"])
    1 <= parameter_index <= length(u_values) || error("parameter index out of range")
    u_value = u_values[parameter_index]
    ground_sector_index = Int(config["observables"]["ground_sector_indices"][parameter_index])
    production_sector = production_sector_specs(torus)[ground_sector_index]
    momentum = production_sector.momentum

    u_milli = round(Int, 1000 * u_value)
    output_path = joinpath(
        abspath(config["output"]["root"]), "ground_chirality",
        "u_$(lpad(u_milli, 5, '0')).jld2",
    )
    if isfile(output_path) && !Bool(config["output"]["overwrite"])
        println("skip existing: ", output_path)
        return
    end

    println("ground chirality U/t=", u_value, ", momentum=", (momentum.nx, momentum.ny))
    block = xdiag_electron_block(torus, 6, 6, momentum)
    hamiltonian = build_hubbard_hamiltonian(torus; t=1.0, U=u_value)
    operator = XDiag.csr_matrix_32(hamiltonian, block)
    solved = trusted_ground(operator, block, config["solver"])

    reference_sector = production_sector_specs(torus)[ground_sector_index]
    reference_id = sector_id(reference_sector.momentum, reference_sector.parity)
    reference_path = joinpath(
        abspath(config["observables"]["reference_root"]), reference_id,
        "u_$(lpad(u_milli, 5, '0')).jld2",
    )
    reference_energy = Float64(load_result(reference_path)["spectrum"]["energy"][1])
    difference = abs(solved.energy - reference_energy)
    difference <= 1e-7 || error("translation-only ground disagrees by $difference")

    timing = @timed measure_chirality_structure(
        solved.vector, block, torus, momentum
    )
    chirality = timing.value
    chi = real(chirality["chirality_mean_per_triangle"])
    projection = momentum == MomentumSector(0, 0) ? [chi] : sort([-chi, chi])

    record = result_record(
        metadata=Dict(
            "name" => "N12_ground_chirality_structure",
            "created_utc" => string(now(UTC)),
            "config_path" => relpath(config_path, pwd()),
        ),
        cluster=Dict("Lx" => 3, "Ly" => 4, "N" => 12, "nup" => 6, "ndn" => 6),
        parameters=Dict("U_over_t" => u_value, "t" => 1.0, "U" => u_value),
        sector=Dict(
            "id" => sector_id(momentum), "nx" => momentum.nx, "ny" => momentum.ny,
            "C2_parity" => 0, "block_dimension" => XDiag.dim(block),
        ),
        solver=Dict(
            "algorithm" => "eigs_lobpcg", "criterion" => solved.solution.criterion,
            "iterations" => solved.solution.niterations,
            "residual_norms" => [Float64(solved.solution.residual_norms[1])],
            "explicit_residual_norms" => [solved.residual],
            "selected_random_seed" => solved.seed, "solve_seconds" => solved.seconds,
        ),
        spectrum=Dict(
            "energy" => [solved.energy], "energy_minus_block_ground" => [0.0],
            "level_index" => [1],
        ),
        observables=Dict(
            "chirality" => chirality,
            "projected_uniform_chirality_eigenvalues" => projection,
            "measurement_seconds" => timing.time,
            "reference_energy_difference" => difference,
        ),
        provenance=Dict(
            "julia_version" => string(VERSION),
            "xdiag_version" => string(pkgversion(XDiag)),
            "slurm_job_id" => get(ENV, "SLURM_JOB_ID", "interactive"),
        ),
    )
    save_result(output_path, record; overwrite=Bool(config["output"]["overwrite"]))
    println("chi mean=", chi, ", measurement seconds=", timing.time)
    println("saved: ", output_path)
end

main(ARGS)
