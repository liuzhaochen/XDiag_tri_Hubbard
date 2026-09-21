using LinearAlgebra
using Printf
using XDiag
using XDiagTriHubbard

function main(args)
    length(args) >= 1 || error(
        "usage: crosscheck_ground_momentum.jl U_OVER_T... [--nup N --ndn N]"
    )

    nup = 6
    ndn = 6
    u_values = Float64[]
    index = 1
    while index <= length(args)
        if args[index] == "--nup"
            index += 1
            nup = parse(Int, args[index])
        elseif args[index] == "--ndn"
            index += 1
            ndn = parse(Int, args[index])
        else
            push!(u_values, parse(Float64, args[index]))
        end
        index += 1
    end
    isempty(u_values) && error("at least one U/t value is required")

    torus = TriangularTorus(3, 4)

    function independent_ground(hamiltonian, momentum)
        block = xdiag_electron_block(torus, nup, ndn, momentum)
        operator = XDiag.csr_matrix_32(hamiltonian, block)
        solution = XDiag.eigs_lobpcg(
            operator,
            block;
            neigs=1,
            guard=2,
            tol=1e-11,
            max_iterations=1000,
            random_seed=1729,
        )
        solution.criterion == "converged" || error(
            "LOBPCG failed at ($(momentum.nx),$(momentum.ny)): $(solution.criterion)"
        )
        vector = XDiag.matrix(solution.eigenvectors; copy=true)[:, 1]
        energy = Float64(solution.eigenvalues[1])
        residual = norm(XDiag.apply(operator, vector) - energy * vector)
        residual <= 1e-8 || error(
            "residual $residual failed at ($(momentum.nx),$(momentum.ny))"
        )
        return (; momentum, energy, residual, dimension=XDiag.dim(block))
    end

    for u_value in u_values
        hamiltonian = build_hubbard_hamiltonian(torus; t=1.0, U=u_value)
        results = [independent_ground(hamiltonian, momentum)
                   for momentum in all_momenta(torus)]
        sort!(results; by=result -> result.energy)
        minimum_energy = results[1].energy
        @printf("U/t=%.8g, Nup=%d, Ndn=%d\n", u_value, nup, ndn)
        for result in results
            @printf(
                "  (%d,%d) E=%.15f E-Emin=%.12g residual=%.3e dim=%d\n",
                result.momentum.nx,
                result.momentum.ny,
                result.energy,
                result.energy - minimum_energy,
                result.residual,
                result.dimension,
            )
        end
        flush(stdout)
    end
end

main(ARGS)
