using LinearAlgebra
using Printf
using XDiag
using XDiagTriHubbard

length(ARGS) == 1 || error("usage: validate_sector_spectra.jl U_OVER_T")
u_value = parse(Float64, ARGS[1])
torus = TriangularTorus(3, 4)
hamiltonian = build_hubbard_hamiltonian(torus; t=1.0, U=u_value)

function low_energies(momentum::MomentumSector; parity=nothing, neigs::Int=6)
    block = xdiag_electron_block(torus, 6, 6, momentum; parity=parity)
    operator = XDiag.csr_matrix_32(hamiltonian, block)
    for attempt in 1:3
        seed = 43 + (attempt - 1) * 104729
        solution = XDiag.eigs_lobpcg(
            operator,
            block;
            neigs=neigs,
            guard=2,
            tol=1e-10,
            max_iterations=500,
            random_seed=seed,
        )
        if length(solution.eigenvalues) >= neigs &&
           length(solution.residual_norms) >= neigs &&
           solution.criterion == "converged" &&
           maximum(solution.residual_norms[1:neigs]) <= 1e-8
            vectors = XDiag.matrix(solution.eigenvectors; copy=true)[:, 1:neigs]
            explicit_residuals = [
                norm(XDiag.apply(operator, vectors[:, level]) -
                     solution.eigenvalues[level] * vectors[:, level])
                for level in 1:neigs
            ]
            orthogonality = maximum(abs, vectors' * vectors - I)
            if maximum(explicit_residuals) <= 1e-8 && orthogonality <= 1e-8
                return Float64.(solution.eigenvalues[1:neigs])
            end
        end
    end
    error("LOBPCG retries exhausted for momentum $(momentum.nx),$(momentum.ny), parity=$parity")
end

partner_pairs = [
    (MomentumSector(0, 1), MomentumSector(0, 3)),
    (MomentumSector(1, 0), MomentumSector(2, 0)),
    (MomentumSector(1, 1), MomentumSector(2, 3)),
    (MomentumSector(1, 2), MomentumSector(2, 2)),
    (MomentumSector(1, 3), MomentumSector(2, 1)),
]

println("partner and parity validation at U/t=", u_value)
for (first_momentum, second_momentum) in partner_pairs
    first_values = low_energies(first_momentum)
    second_values = low_energies(second_momentum)
    difference = maximum(abs, first_values - second_values)
    @printf("  (%d,%d) <-> (%d,%d): max |dE| = %.3e\n",
            first_momentum.nx, first_momentum.ny,
            second_momentum.nx, second_momentum.ny, difference)
    difference <= 1e-8 || error("momentum-partner spectrum mismatch")
    flush(stdout)
end

for momentum in (MomentumSector(0, 0), MomentumSector(0, 2))
    translation_values = low_energies(momentum)
    even_values = low_energies(momentum; parity=1)
    odd_values = low_energies(momentum; parity=-1)
    reconstructed = sort(vcat(even_values, odd_values))[1:length(translation_values)]
    difference = maximum(abs, translation_values - reconstructed)
    @printf("  (%d,%d) translation <-> A/B irrep union: max |dE| = %.3e\n",
            momentum.nx, momentum.ny, difference)
    difference <= 1e-8 || error("C2-resolved spectrum mismatch")
    flush(stdout)
end

println("validation passed at U/t=", u_value)
