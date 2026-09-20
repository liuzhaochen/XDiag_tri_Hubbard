using LinearAlgebra
using Printf
using XDiag
using XDiagTriHubbard

torus = TriangularTorus(3, 4)
momentum = MomentumSector(1, 0)
block = xdiag_electron_block(torus, 6, 6, momentum)
hamiltonian = build_hubbard_hamiltonian(torus; t=1.0, U=8.0)

println("benchmark sector: (1,0), dimension=", XDiag.dim(block), ", U/t=8")
flush(stdout)

csr_timing = @timed XDiag.csr_matrix_32(hamiltonian, block)
operator = csr_timing.value
println("CSR32: nnz=", length(operator.data), ", seconds=", csr_timing.time,
        ", bytes=", Base.summarysize(operator))
flush(stdout)

function diagnostics(operator, solution, requested::Int)
    vectors = XDiag.matrix(solution.eigenvectors; copy=true)[:, 1:requested]
    values = solution.eigenvalues[1:requested]
    residuals = [norm(XDiag.apply(operator, vectors[:, level]) -
                      values[level] * vectors[:, level])
                 for level in 1:requested]
    orthogonality = maximum(abs, vectors' * vectors - I)
    return values, residuals, orthogonality
end

requested = 6
lanczos_timing = @timed XDiag.eigs_lanczos(
    operator,
    block;
    neigvals=requested,
    precision=1e-10,
    max_iterations=500,
    deflation_tol=1e-8,
    random_seed=43,
)
lanczos = lanczos_timing.value
lanczos_values, lanczos_residuals, lanczos_orthogonality =
    diagnostics(operator, lanczos, requested)
println("Lanczos: criterion=", lanczos.criterion,
        ", iterations=", lanczos.niterations,
        ", seconds=", lanczos_timing.time)
println("Lanczos energies: ", lanczos_values)
println("Lanczos residuals: ", lanczos_residuals)
println("Lanczos orthogonality error: ", lanczos_orthogonality)
flush(stdout)

lobpcg_timing = @timed XDiag.eigs_lobpcg(
    operator,
    block;
    neigs=requested,
    guard=2,
    tol=1e-10,
    max_iterations=500,
    random_seed=43,
)
lobpcg = lobpcg_timing.value
lobpcg_values, explicit_lobpcg_residuals, lobpcg_orthogonality =
    diagnostics(operator, lobpcg, requested)
println("LOBPCG: criterion=", lobpcg.criterion,
        ", iterations=", lobpcg.niterations,
        ", seconds=", lobpcg_timing.time)
println("LOBPCG energies: ", lobpcg_values)
println("LOBPCG reported residuals: ", lobpcg.residual_norms[1:requested])
println("LOBPCG explicit residuals: ", explicit_lobpcg_residuals)
println("LOBPCG orthogonality error: ", lobpcg_orthogonality)
println("maximum Lanczos/LOBPCG energy difference: ",
        maximum(abs, lanczos_values - lobpcg_values))

lanczos_rejected = maximum(lanczos_residuals) > 1e-8 ||
                   lanczos_orthogonality > 1e-8
lanczos_rejected || error("benchmark was expected to expose unsafe multilevel Lanczos")
abs(lanczos_values[1] - lobpcg_values[1]) <= 1e-10 ||
    error("trusted ground values disagree")
lobpcg.criterion == "converged" || error("LOBPCG did not converge")
maximum(lobpcg.residual_norms[1:requested]) <= 1e-8 ||
    error("LOBPCG residual gate failed")
lobpcg_orthogonality <= 1e-8 || error("LOBPCG orthogonality gate failed")
println("decision: reject multilevel Lanczos; accept guarded LOBPCG")
