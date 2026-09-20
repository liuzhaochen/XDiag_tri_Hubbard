using Printf
using XDiag
using XDiagTriHubbard

torus = TriangularTorus(3, 4)
println("Julia version: ", VERSION)
println("XDiag version: ", pkgversion(XDiag))
println("raw fixed-particle dimension: ", binomial(12, 6)^2)
println("production blocks:")

for sector in production_sector_specs(torus)
    momentum = sector.momentum
    block = xdiag_electron_block(
        torus, 6, 6, momentum; parity=sector.parity
    )
    parity_label = sector.parity === nothing ? "-" : string(sector.parity)
    @printf("  (%d,%d) C2=%2s star=%d dimension=%d\n",
            momentum.nx, momentum.ny, parity_label, length(sector.star), XDiag.dim(block))
end
