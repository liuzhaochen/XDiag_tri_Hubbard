using Dates
using LinearAlgebra
using TOML
using XDiag
using XDiagTriHubbard

function main(args)
length(args) == 3 || error(
    "usage: run_scan.jl CONFIG SECTOR_INDEX PARAMETER_INDEX"
)
config_path = abspath(args[1])
sector_index = parse(Int, args[2])
parameter_index = parse(Int, args[3])
config = TOML.parsefile(config_path)

cluster_config = config["cluster"]
torus = TriangularTorus(Int(cluster_config["Lx"]), Int(cluster_config["Ly"]))
length(torus) == Int(cluster_config["N"]) || error("cluster size mismatch")
nup, ndn = Int(cluster_config["nup"]), Int(cluster_config["ndn"])
nup == 6 && ndn == 6 || error("production requires Nup=Ndn=6")

all_sectors = production_sector_specs(torus)
selected_indices = haskey(config["sectors"], "indices") ?
    Int.(config["sectors"]["indices"]) : collect(eachindex(all_sectors))
all(index -> 1 <= index <= length(all_sectors), selected_indices) ||
    error("selected sector index out of range")
length(Set(selected_indices)) == length(selected_indices) ||
    error("selected sector indices contain duplicates")
sectors = all_sectors[selected_indices]
length(sectors) == Int(config["sectors"]["expected_count"]) ||
    error("selected sector count mismatch")
1 <= sector_index <= length(sectors) || error("sector index out of range")
sector = sectors[sector_index]
momentum, parity = sector.momentum, sector.parity
identifier = sector_id(momentum, parity)

u_values = Float64.(config["trajectory"]["values"])
1 <= parameter_index <= length(u_values) || error("parameter index out of range")
u_value = u_values[parameter_index]
t_value = Float64(config["trajectory"]["t"])

output_root = abspath(config["output"]["root"])
u_milli = round(Int, 1000 * u_value)
output_path = joinpath(
    output_root, identifier, "u_$(lpad(u_milli, 5, '0')).jld2"
)
overwrite = Bool(config["output"]["overwrite"])
if isfile(output_path) && !overwrite
    println("skip existing: ", output_path)
    exit()
end

println("sector=", identifier, ", U/t=", u_value)
block_timing = @timed xdiag_electron_block(
    torus, nup, ndn, momentum; parity=parity
)
block = block_timing.value
block_dimension = XDiag.dim(block)
println("block dimension=", block_dimension,
        ", construction seconds=", block_timing.time)

hamiltonian = build_hubbard_hamiltonian(torus; t=t_value, U=u_value * t_value)
solver_config = config["solver"]
solver_config["algorithm"] == "eigs_lobpcg" || error("unsupported solver")
solver_config["operator_mode"] == "csr32" || error("unsupported operator mode")

operator_timing = @timed XDiag.csr_matrix_32(hamiltonian, block)
operator = operator_timing.value
println("CSR32 nnz=", length(operator.data),
        ", construction seconds=", operator_timing.time)
flush(stdout)

requested = Int(solver_config["neigenvalues"])
residual_limit = Float64(solver_config["residual_acceptance_limit"])
orthogonality_limit = Float64(solver_config["orthogonality_tolerance"])
base_seed = Int(solver_config["random_seed"])
retry_count = Int(solver_config["retry_count"])
seed_stride = Int(solver_config["retry_seed_stride"])

solution = nothing
solve_seconds = 0.0
selected_seed = 0
failed_attempts = String[]
explicit_residuals = Float64[]
orthogonality_error = Inf
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
    if length(candidate.eigenvalues) >= requested &&
       length(candidate.residual_norms) >= requested
        vectors = XDiag.matrix(candidate.eigenvectors; copy=true)[:, 1:requested]
        values = candidate.eigenvalues[1:requested]
        explicit_residuals = [
            norm(XDiag.apply(operator, vectors[:, level]) -
                 values[level] * vectors[:, level]) for level in 1:requested
        ]
        orthogonality_error = maximum(abs, vectors' * vectors - I)
        residual_ok = maximum(candidate.residual_norms[1:requested]) <= residual_limit &&
                      maximum(explicit_residuals) <= residual_limit
        orthogonality_ok = orthogonality_error <= orthogonality_limit
        if candidate.criterion == "converged" && residual_ok && orthogonality_ok
            solution = candidate
            selected_seed = seed
            break
        end
        push!(failed_attempts,
              "seed=$seed criterion=$(candidate.criterion) " *
              "reported=$(maximum(candidate.residual_norms[1:requested])) " *
              "explicit=$(maximum(explicit_residuals)) orthogonality=$orthogonality_error")
    else
        push!(failed_attempts, "seed=$seed returned too few states")
    end
end
solution === nothing && error("LOBPCG retries exhausted: $(join(failed_attempts, "; "))")

energies = Float64.(solution.eigenvalues[1:requested])
residuals = Float64.(solution.residual_norms[1:requested])
println("energies=", energies)
println("maximum residual=", maximum(residuals),
        ", orthogonality error=", orthogonality_error,
        ", solve seconds=", solve_seconds)

function git_revision()
    try
        return readchomp(pipeline(`git rev-parse HEAD`; stderr=devnull))
    catch
        return "unavailable"
    end
end

record = result_record(
    metadata=Dict(
        "name" => config["name"],
        "created_utc" => string(now(UTC)),
        "config_path" => relpath(config_path, pwd()),
    ),
    cluster=Dict(
        "Lx" => torus.Lx,
        "Ly" => torus.Ly,
        "N" => length(torus),
        "nup" => nup,
        "ndn" => ndn,
        "bond_count" => length(nearest_neighbor_bonds(torus)),
    ),
    parameters=Dict(
        "U_over_t" => u_value,
        "t" => t_value,
        "U" => u_value * t_value,
    ),
    sector=Dict(
        "id" => identifier,
        "nx" => momentum.nx,
        "ny" => momentum.ny,
        "Ky" => 2pi * momentum.ny / torus.Ly,
        "C2_parity" => parity === nothing ? 0 : parity,
        "star_nx" => [k.nx for k in sector.star],
        "star_ny" => [k.ny for k in sector.star],
        "star_size" => length(sector.star),
        "block_dimension" => block_dimension,
    ),
    solver=Dict(
        "algorithm" => "eigs_lobpcg",
        "operator_mode" => "csr32",
        "criterion" => solution.criterion,
        "iterations" => solution.niterations,
        "requested_eigenvalues" => requested,
        "guard_eigenvalues" => Int(solver_config["guard_eigenvalues"]),
        "precision" => Float64(solver_config["precision"]),
        "residual_acceptance_limit" => residual_limit,
        "residual_norms" => residuals,
        "explicit_residual_norms" => explicit_residuals,
        "orthogonality_error" => orthogonality_error,
        "selected_random_seed" => selected_seed,
        "failed_attempts" => failed_attempts,
        "block_construction_seconds" => block_timing.time,
        "operator_construction_seconds" => operator_timing.time,
        "solve_seconds" => solve_seconds,
        "operator_nonzeros" => length(operator.data),
    ),
    spectrum=Dict(
        "energy" => energies,
        "energy_minus_block_ground" => energies .- minimum(energies),
        "level_index" => collect(1:requested),
    ),
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
