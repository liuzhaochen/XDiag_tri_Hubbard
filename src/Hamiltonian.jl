import XDiag

"""Construct the direct triangular-lattice Hubbard Hamiltonian.

XDiag's `Hop` is already Hermitian, includes both spin species, and implements the
minus sign in `-t(cdag_i*c_j+h.c.)`. Add it once per undirected bond.
"""
function build_hubbard_hamiltonian(
    torus::TriangularTorus;
    t::Real=1.0,
    U::Real=0.0,
)
    hamiltonian = XDiag.OpSum()
    for (i, j) in nearest_neighbor_bonds(torus)
        hamiltonian += Float64(t) * XDiag.Op("Hop", [i, j])
    end
    hamiltonian += Float64(U) * XDiag.Op("HubbardU")
    return hamiltonian
end
