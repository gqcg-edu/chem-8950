## Restricted Hartree-Fock Theory

In the previous set of notes, we took our energy expression 
\\[E =  \sum\limits_{i}^N \langle\psi_i^i|\hat{h}(i)|\psi_i^i \rangle + \sum\limits_{i<j}^N  \langle \psi_i^i\psi_j^j|\hat{g}(i,j)|\psi_i^i\psi_j^j \rangle - \langle \psi_i^i\psi_j^j|\hat{g}(i,j)|\psi_j^i\psi_i^j \rangle \\]

and derived the closed-shell special case, where we have integrated out all our spin functions and now just have an expression in terms of just spatial orbitals \$\phi\$:


\\[E = 2 \sum\limits_{i}^{MO} \langle \phi_i | \hat{h}(i) | \phi_i \rangle +  \sum\limits_{i}^{MO} \sum\limits_{j}^{MO} 2 \langle \phi_i \phi_j | \hat{g}(i,j) | \phi_i \phi_j \rangle - \langle \phi_i \phi_j | \hat{g}(i,j) | \phi_j \phi_i \rangle\\]

We have not specified what our oribtals are. In principle, they could be any function of 3 spatial variables. Further, it is not immediately obvious whether a given energy from one chosen set of orbitals is better from _another_ energy given by _another_ chosen set of orbitals. Well, the Variational Principle states that any chosen trial wavefunction (and thus, set of orbitals) will always overestimate (be "bounded below") by the true ground state energy. So, varying the orbitals such that the energy is lowest will give the best wavefunction, and best energy approximation. So, some change \$\delta\$ (ideally a decrease) in the energy can be induced by some small change \$\delta\$ in the orbitals, which we might unrigorously express as:

\\[ \delta E = 2 \sum\limits_{i}^{MO} \langle \delta \phi_i | \hat{h}(i) | \delta \phi_i \rangle +  \sum\limits_{i}^{MO} \sum\limits_{j}^{MO} 2 \langle \delta \phi_i \delta\phi_j | \hat{g}(i,j) | \delta\phi_i \delta\phi_j \rangle - \langle \delta\phi_i \delta\phi_j | \hat{g}(i,j) |\delta \phi_j \delta\phi_i \rangle \\]

where it is presumed that different orbitals will have different changes \$\delta\$. But, we cannot just vary the orbitals however we wish. This energy expression was derived under the assumption that the orbitals are orthonormal, so we need to keep them that way. Thus, we have a **minimization** problem of the energy under some **constraint** that the orbitals are kept orthonormal.

The equations which can be used satisfy these conditions are known as the Hartree-Fock equations. We will not derive the Hartree-Fock equations here. I recommend the following resources:
    1. Roothaan's 1951 paper, "New Developments in Molcular Orbital Theory"
    2. Dr. Andreas Copan's notes from the 2017 edition of this course (on GitHub)
    3. Szabo and Ostlund, "Modern Quantum Chemistry"

All of these derivations would almost certainly require you to brush-up on the calculus of variations, in particular the method of Lagrange multipliers. I would only recommend doing this if you find it likely you will use these routines in future work. 

### Hartree-Fock Equations
First we need to introduce the "Coulomb" \\(\hat{J}_j\\) and "Exchange" \$\hat{K}_j\$ *operators*. These are really only defined for aesthetic purposes; they make our equations simpler to express. We define them in terms of how they act on an orbital:

\\[ \hat{J}_j |\phi^\mu \rangle = \langle \phi_j^\nu |\hat{g}(\mu,\nu)| \phi_j^\nu \phi^\mu \rangle  \\]

\\[ \hat{K}_j |\phi^\mu \rangle = \langle \phi_j^\nu |\hat{g}(\mu,\nu)| \phi^\nu \phi_j^\mu \rangle  \\]

The exchange operator is without a doubt a very weird operator; it rips the electronic coordinates out of the orbital its operating on, and stuffs in some other electron's coordinates instead. Under this new notation, our energy expression looks like the following:

\\[E = 2 \sum\limits_{i}^{MO} \langle \phi_i | \hat{h}(i) | \phi_i \rangle +  \sum\limits_{i}^{MO} \sum\limits_{j}^{MO} 2 \langle \phi_i  | \hat{J}_j | \phi_i \rangle - \langle \phi_i | \hat{K}_j | \phi_i \rangle \\]

Now, every operator is a "one-electron operator", and by this we mean that each operator gives a piece of the energy when we take the expectation value over some one-electron wavefunction \$\phi_i\$.

In the Hartree-Fock equation derivation, one finds that the best set of molecular orbitals satisfy the following equations:

\\[ \left[\hat{h} + \sum\limits_j (2 \hat{J}_j - \hat{K}_j) \right] \phi_i = \epsilon_i \phi_i  \\]

The term in square brackets is referred to as the Fock Operator \$\hat{F}\$, so we just write

\\[\hat{F} \phi_i = \epsilon_i \phi_i \\]

These are the Hartree-Fock (HF) equations. This simple set of equations for a set of orbitals \\(\{ \phi_i \}\\) represent a great deal of buried complexity. The MO's which give the best Slater determinant (and lowest energy) are all eigenfunctions of the hermitian operator \$\hat{F}\$, which itself is defined _in terms of these MO's_ (look at the equations of our coulomb and exchange operators).

Thus, solving these equations requires trial and error. Pick a set of \$\phi_i\$, build \$\hat{F}\$, solve the equations for the \$n\$ lowest eigenvalues, compare the newly obtained \$\phi_i\$'s to the old ones. One repeats this procedure until the \$\phi_i\$'s converge. Hence, this procedure is commonly called the _self-consisten field_ procedure. Our 'field' of energy interactions described by the Fock operator must stabilize to the optimal solution when the orbitals are improved. **The self-consistent orbitals are the best MO's from which to construct a Slater determinant wavefunction. That Slater determinant wavefunction is the best possible single Slater determinant wavefunction obtainable. The corresponding energy is the best possible approximate ground state energy, and it is bounded below by the true electronic energy of the system.**

The HF equations are quite expensive to solve. For atoms, it is a bit easier due to the spherical symmetry, and since we kind of know what the one-electron wavefunctions should look like. For molecules, solving the HF equations is pretty much impossible. **As a result, we have to use an approximation to the best MO's, since obtaining them analytically by the HF method is unfeasible in the fast majority of cases.**


## The Roothaan-Hall Equations 

To apply the Hartree-Fock method to molecules, we approximate the best MO's which solve the Hartree-Fock equations by a _linear combination of atomic orbitals_ (LCAO), which form our MO's. We will call these LCAO-MO's, to distinguish them from the MO's of Hartree-Fock theory.

\\[\phi_i = \sum\limits_p \chi_p C_{pi} \\]

