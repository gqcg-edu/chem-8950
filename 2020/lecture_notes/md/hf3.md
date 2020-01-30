## Restricted Hartree-Fock Theory

In the previous set of notes, we took our energy expression 
\\[E =  \sum\limits_{i}^N \langle\psi_i^i|\hat{h}(i)|\psi_i^i \rangle + \sum\limits_{i<j}^N  \langle \psi_i^i\psi_j^j|\hat{g}(i,j)|\psi_i^i\psi_j^j \rangle - \langle \psi_i^i\psi_j^j|\hat{g}(i,j)|\psi_j^i\psi_i^j \rangle \\]

and derived the closed-shell special case, where we have integrated out all our spin functions and now just have an expression in terms of just spatial orbitals \\(\phi\\) and sums over \\(N/2\\) electrons (the number of doubly-occupied MO's):


\\[E = 2 \sum\limits_{i}^{N/2} \langle \phi_i | \hat{h}(i) | \phi_i \rangle +  \sum\limits_{i}^{N/2} \sum\limits_{j}^{N/2} 2 \langle \phi_i \phi_j | \hat{g}(i,j) | \phi_i \phi_j \rangle - \langle \phi_i \phi_j | \hat{g}(i,j) | \phi_j \phi_i \rangle\\]

We have not specified what our oribtals are. In principle, they could be any function of 3 spatial variables. Further, it is not immediately obvious whether a given energy from one chosen set of orbitals is better than _another_ energy given by _another_ chosen set of orbitals. Well, the Variational Principle states that any chosen trial wavefunction (and in our case, a set of orbitals) will always overestimate (be "bounded below") by the true ground state energy. So, varying the orbitals such that the energy is lowest will give the best wavefunction, and best energy approximation. Some change \\(\delta\\) (ideally a decrease) in the energy can be induced by some small change \\(\delta\\) in the orbitals, which we might unrigorously express as:

\\[ \delta E = 2 \sum\limits_{i}^{N/2} \langle \delta \phi_i | \hat{h}(i) | \delta \phi_i \rangle +  \sum\limits_{i}^{N/2} \sum\limits_{j}^{N/2} 2 \langle \delta \phi_i \delta\phi_j | \hat{g}(i,j) | \delta\phi_i \delta\phi_j \rangle - \langle \delta\phi_i \delta\phi_j | \hat{g}(i,j) |\delta \phi_j \delta\phi_i \rangle \\]

where it is presumed that different orbitals will have different changes \\(\delta\\). But, we cannot just vary the orbitals however we wish. This energy expression was derived under the assumption that the orbitals are orthonormal, so we need to keep them that way. Thus, we have a **minimization** problem of the energy under some **constraint** that the orbitals are kept orthonormal.

The equations which can be used satisfy these conditions are known as the Hartree-Fock equations. We will not derive the Hartree-Fock equations here. I recommend the following resources:
    1. Roothaan's 1951 paper, "New Developments in Molcular Orbital Theory"
    2. Dr. Andreas Copan's notes from the 2017 edition of this course (on GitHub)
    3. Szabo and Ostlund, "Modern Quantum Chemistry"

All of these derivations would almost certainly require you to brush-up on the calculus of variations, in particular the method of Lagrange multipliers. I would only recommend doing this if you find it likely you will use these routines in future work. 

### Hartree-Fock Equations
First we need to introduce the "coulomb" \\(\hat{J}_j\\) and "exchange" \\(\hat{K}_j\\) *operators*. These are really only defined for aesthetic purposes; they make our equations simpler to express. We define them in terms of how they act on an orbital:

\\[ \hat{J}_j |\phi^\mu \rangle = \langle \phi_j^\nu |\hat{g}(\mu,\nu)| \phi_j^\nu \phi^\mu \rangle  \\]

\\[ \hat{K}_j |\phi^\mu \rangle = \langle \phi_j^\nu |\hat{g}(\mu,\nu)| \phi^\nu \phi_j^\mu \rangle  \\]

The exchange operator is without a doubt a very weird operator; it rips the electronic coordinates out of the orbital its operating on, and stuffs in some other electron's coordinates instead. Under this new notation, our energy expression looks like the following:

\\[E = 2 \sum\limits_{i}^{N/2} \langle \phi_i | \hat{h}(i) | \phi_i \rangle +  \sum\limits_{i}^{N/2} \sum\limits_{j}^{N/2} 2 \langle \phi_i  | \hat{J}_j | \phi_i \rangle - \langle \phi_i | \hat{K}_j | \phi_i \rangle \\]

Now, every operator is a "one-electron operator", and by this we mean that each operator gives a piece of the energy when we take the expectation value over some one-electron wavefunction \\(\phi_i\\).

In the Hartree-Fock equation derivation, one finds that the best set of molecular orbitals satisfy the following equations:

\\[ \left[\hat{h} + \sum\limits_j (2 \hat{J}_j - \hat{K}_j) \right] \phi_i = \epsilon_i \phi_i  \\]

The term in square brackets is referred to as the Fock Operator \\(\hat{f}\\), so we just write

\\[\hat{f} \phi_i = \epsilon_i \phi_i \\]

These are the Hartree-Fock (HF) equations. This simple set of equations for a set of orbitals \\(\{ \phi_i \}\\) represent a great deal of buried complexity. The MO's which give the best Slater determinant (and lowest energy) are all eigenfunctions of the hermitian operator \\(\hat{f}\\), which itself is defined _in terms of these MO's_ (look at the equations of our coulomb and exchange operators).

Thus, solving these equations requires trial and error. Pick a set of \\(\phi_i\\), build \\(\hat{f}\\), solve the equations for the \\(n\\) lowest eigenvalues (orbital energies \\(\epsilon_i\\)), compare the newly obtained \\(\phi_i\\)'s to the old ones. One repeats this procedure until the \\(\phi_i\\)'s converge. Hence, this procedure is commonly called the _self-consistent field_ procedure. Our 'field' of energy interactions described by the Fock operator must stabilize to the optimal solution as the orbitals are improved. **The self-consistent (converged) orbitals are the best MO's from which to construct a Slater determinant wavefunction. That Slater determinant wavefunction is the best possible single Slater determinant wavefunction obtainable. The corresponding energy is the best possible approximate ground state energy, and it is bounded below by the true electronic energy of the system.**

The HF equations are quite expensive to solve. For atoms, it is a bit easier due to the spherical symmetry, and since we kind of know what the one-electron wavefunctions should look like. For molecules, solving the HF equations is pretty much impossible. **As a result, we have to use an approximation to the best MO's, since obtaining them analytically by the HF method is unfeasible in the fast majority of cases.**


## The Roothaan-Hall Equations 

To apply the Hartree-Fock method to molecules, we approximate the best MO's which solve the Hartree-Fock equations by a _linear combination of atomic orbitals_ (LCAO), which form our MO's. We will denote them with the verbose acronym "LCAO-MO", to distinguish them from the 'true' MO's of Hartree-Fock theory. One can think of the AO's (\\(\chi\\)) as a _basis_ for representing the LCAO-MO's. In the language of our state vector discussions: a one electron wavefunction (a LCAO-MO, \\(\phi\\)) can be represented by a linear combination of _base states_ (our AO basis functions) with certain coefficients. Mathematically,

\\[\phi_i = \sum\limits_q \chi_q C_{qi} \\]


We can approximate the Hartree-Fock MO's with the LCAO-MO's by inserting the expansion above into the HF equations:

\\[\hat{f} \phi_i = \epsilon_i \phi_i \\]

\\[\hat{f} \sum\limits_q \chi_q C_{qi} = \epsilon_i \sum\limits_q \chi_q C_{qi} \\]

Left-multiplying by some other AO \\(\chi_p^*\\) and integrating on both sides gives

\\[ \sum\limits_q \langle \chi_p \mid \hat{f}\mid \chi_q \rangle C_{qi} = \epsilon_i \sum\limits_q \langle \chi_p \mid \chi_q \rangle C_{qi} \\]

Defining the integrals above \\(\langle \chi_p \mid \hat{f}\mid \chi_q \rangle \\) and \\(\langle \chi_p \mid \chi_q \rangle\\) as \\(F_{pq} \\) and \\(S_{pq}\\),

\\[ \sum\limits_q F_{pq} C_{qi} = \epsilon_i \sum\limits_q S_{pq} C_{qi} \\]


The above equations are valid for a particular \\(\phi_i\\). If we want a single equation which determines **all** \\(\phi_i\\)'s in terms of our set of \\(\chi \\)'s, we can reformulate the above expression in matrix notation
\\[\boldsymbol{F}_{pq}\boldsymbol{C}_{qi} = \boldsymbol{S}_{pq}\boldsymbol{C}_{qi}\boldsymbol{\epsilon}_{i}    \\]
so that for \\(m\\) AO basis functions, \\(\boldsymbol{C}_{qi} \\) is an \\(m \times m\\) matrix with each column vector containing the expansion coefficients for a LCAO-MO \\(\phi_i\\). The Fock matrix \\(\boldsymbol{F}_{pq} \\) and overlap matrix \\(\boldsymbol{S}_{pq} \\) are defined as
\begin{align}
\boldsymbol{F}_{pq} =
\begin{pmatrix}
\langle \chi_1 \mid \hat{f} \mid \chi_1 \rangle & \dots & \langle \chi_1 \mid \hat{f} \mid \chi_m \rangle \\ 
\vdots & \ddots & \vdots \\
\langle \chi_m \mid \hat{f} \mid \chi_1 \rangle & \dots & \langle \chi_m \mid \hat{f} \mid \chi_m \rangle \\ 
\end{pmatrix}
\end{align}

\begin{align}
\boldsymbol{S}_{pq} =
\begin{pmatrix}
\langle \chi_1 \mid \chi_1 \rangle & \dots & \langle \chi_1  \mid \chi_m \rangle \\ 
\vdots & \ddots & \vdots \\
\langle \chi_m  \mid \chi_1 \rangle & \dots & \langle \chi_m   \mid \chi_m \rangle \\ 
\end{pmatrix}
\end{align}

Having defined these matrices explicitly, we will drop the subscripts from here on, and just write
\\[\boldsymbol{F}\boldsymbol{C} = \boldsymbol{S}\boldsymbol{C}\boldsymbol{\epsilon}  \\]
The above is known as the Roothaan-Hall equations. Clemens Roothaan and George Hall independently derived and published this scheme in 1951. Not to play favorites or anything, but Roothaan's paper is far more intelligible and elegant than Hall's, so if you're the kind of person who likes to learn from the source, read Roothaan's paper.

### Fock Matrix Elements
It is worth investigating what exactly the \\(\langle \chi_p \mid \hat{f}\mid \chi_q \rangle\\) integrals obtained above are. Using our definition of the Fock operator we obtain:

\\[\langle \chi_p \mid \hat{f}\mid \chi_q \rangle = \langle \chi_p \mid \hat{h} + \sum\limits_j^{N/2} 2 \hat{J}_j - \hat{K}_j \mid \chi_q \rangle  \\]
where the summation is over occupied LCAO-MO's, the number of electrons \\(N\\) divided by two. We can separate each operator into its own bra-ket:

\\[\langle \chi_p \mid \hat{f}\mid \chi_q \rangle = \langle \chi_p \mid \hat{h} \mid \chi_q \rangle + \sum\limits_j^{N/2} 2 \langle \chi_p \mid \hat{J}_j \mid \chi_q \rangle -  \langle \chi_p \mid \hat{K}_j \mid \chi_q \rangle  \\]

Applying the definition of the couloumb and exchange operators we obtain

\\[\langle \chi_p \mid \hat{f}\mid \chi_q \rangle = \langle \chi_p \mid \hat{h} \mid \chi_q \rangle + \sum\limits_j^{N/2} 2 \langle \chi_p \phi_j \mid \hat{g} \mid \chi_q \phi_j \rangle -  \langle \chi_p \phi_j \mid \hat{g} \mid \phi_j \chi_q \rangle  \\]
A brief reminder: we are not being explicit here about which electron's coordinates go where in the two electron integrals above. We are using the implied order of \\(\langle 1,2 \mid \hat{g} \mid 1,2 \rangle \\) for the two electrons. I point this out because if one ignores all notion of electron coordinates, the coulomb and exchange parts above don't look mathematically distinct, but they are. Oftentimes people represent this by writing things like \\(\langle \chi_p(1) \phi_j(2) | \hat{g}(1,2) | \chi_q(1) \phi_j(2) \rangle\\) and \\(\langle \chi_p(1) \phi_j(2) | \hat{g}(1,2) | \phi_j(1) \chi_p(2) \rangle\\) to remind people that the functions have different coordinates, but this makes things very hard to read. Anyways...

We can now expand each \\( \phi_j\\) in the AO basis like we did before. Here it is proper to distinguish between the \\( \phi_j\\) in a bra and \\( \phi_j\\) in a ket by doing two distinct expansions, since our orbitals are not necessarily real (though they usually are):

\\[\langle \phi_j \mid = \sum\limits_r \langle \chi_r \mid C_{rj}^* \quad \quad \mid \phi_j \rangle = \sum\limits_s \mid \chi_s \rangle C_{sj} \\]

Plugging these into our Fock matrix element equation above we obtain:

\\[\langle \chi_p \mid \hat{f}\mid \chi_q \rangle= \langle \chi_p \mid \hat{h} \mid \chi_q \rangle + \sum\limits_j^{N/2}\sum\limits_r^m \sum\limits_s^m C_{rj}^* C_{sj} \left[ 2 \langle \chi_p \chi_r \mid \hat{g} \mid \chi_q \chi_s \rangle - \langle \chi_p \chi_r \mid \hat{g} \mid \chi_s \chi_q \rangle \right]\\]

Defining the density matrix as \\(D_{rs} = \sum\limits_j^{N/2} C_{rj}^* C_{sj} \\), we have:

\\[\langle \chi_p \mid \hat{f}\mid \chi_q \rangle= \langle \chi_p \mid \hat{h} \mid \chi_q \rangle + \sum\limits_r^m \sum\limits_s^m D_{rs} \left[ 2 \langle \chi_p \chi_r \mid \hat{g} \mid \chi_q \chi_s \rangle - \langle \chi_p \chi_r \mid \hat{g} \mid \chi_s \chi_q \rangle \right]\\]

Remembering that our one electron operator \\(\hat{h}\\) is a sum of electron kinetic energy (\\(  \hat{T} = \frac{1}{2} \nabla^2 \\)) and electron-nucleus attraction operators  (\\(\hat{V} = \sum\limits_A \frac{Z_A}{\mid\boldsymbol{r} - \boldsymbol{R}_A \mid}  \\)) we might write the above as the following:

\\[\langle \chi_p \mid \hat{f}\mid \chi_q \rangle= \langle \chi_p \mid \hat{T} \mid \chi_q \rangle + \langle \chi_p \mid \hat{V} \mid \chi_q \rangle + \sum\limits_r^m \sum\limits_s^m D_{rs} \left[ 2 \langle \chi_p \chi_r \mid \hat{g} \mid \chi_q \chi_s \rangle - \langle \chi_p \chi_r \mid \hat{g} \mid \chi_s \chi_q \rangle \right]\\]

If you have coded RHF using Psi4 and NumPy, this should hopefully look familiar. A single Fock matrix element is just a single kinetic energy integral over two AO basis functions, a single potential energy integral over two AO basis functions, and a sum over weighted (by our MO coefficients; our density matrix elements) coulomb and exchange contributions from our set of two-electron integrals. The whole Fock matrix can be built from the matrices which collect all possible kinetic and potential integrals over our AO basis functions, and our two-electron integrals 'contracted' with our density matrix. Keep in mind when I say 'integral', I mean the scalar-value result **of** the integral. Every bra-ket term above evaluates to a scalar, floating point number.

In summary, each element of the Fock matrix is composed of a simple sum of integrals over AO basis functions. The two-electron integral terms are each contracted with the LCAO-MO coefficients. If you have some basis of AO's (a set of \\(\chi\\)'s, known as _basis functions_) for which the above integrals are solvable, and some set of coefficients (perhaps some smart guess?), you can build the Fock matrix. Note this would require computing every integral over every combination of basis functions. Thus, for \\(m\\) basis functions, the one-electron integrals would be collected in an \\(m \times m\\) matrix and the two-electron integrals would be collected in a \\(m \times m \times m \times m\\) array, a rank-4 tensor. From these matrices and two-electron integral array, you can build the Fock matrix, and can then proceed to solve the Roothaan-Hall equations iteratively. 

### Solving the Roothaan-Hall equations

The equation 
\\[\boldsymbol{F}\boldsymbol{C} = \boldsymbol{S}\boldsymbol{C}\boldsymbol{\epsilon}  \\]
is a standard eigenvalue problem if \\(\boldsymbol{S} \\) is an identity matrix. If this is the case, the Roothaan-Hall equations are easy to solve; just diagonalize \\(\boldsymbol{F}\\)! Also, if \\(\boldsymbol{S} \\) is an identity matrix, our basis functions are orthonormal. Recall that our energy expression derivation and the Hartree-Fock equations derivation depended on the basis functions being orthonormal. So 'orthogonalizing' our basis is not only _sufficient_ to make the Roothaan-Hall equations easy to solve, it's also _necessary_ because our equation derivations depended on orbital orthonormality. So, we seek some unitary matrix \\(\boldsymbol{U}\\) such that the transformation performed by this matrix changes our basis to an orthonormal one:

\\[\boldsymbol{U}^T\boldsymbol{S}\boldsymbol{U} = \boldsymbol{I} \\]

There's two such transformations \\(\boldsymbol{U} \\) that are often used. The first is known as the **canonical orthogonalization**, which is found by diagonalizing the overlap matrix:

\\[\boldsymbol{V}^T\boldsymbol{S}\boldsymbol{V} = \boldsymbol{\Lambda} \\]

\\[\boldsymbol{U} = \boldsymbol{V} \boldsymbol{\Lambda}^{-1/2} \\]

\\[\boldsymbol{U}^T \boldsymbol{S} \boldsymbol{U} = \boldsymbol{\Lambda}^{-1/2} \boldsymbol{V}^T \boldsymbol{S} \boldsymbol{V} \boldsymbol{\Lambda}^{-1/2} = \boldsymbol{\Lambda}^{-1/2} \boldsymbol{\Lambda} \boldsymbol{\Lambda}^{-1/2} = \boldsymbol{I}  \\]
The canonical orthogonalization carries the advantage that we can identify linear dependency in the basis set (an eigenvalue in \\(\boldsymbol{\Lambda}\\) would be near zero), and **remove** the linear dependency by just removing the corresponding column of \\( \boldsymbol{V}\\) by setting its eigenvalue partner to zero.

The second \\(\boldsymbol{U}\\) used often to orthogonalize our basis is known as the **symmetric orthogonalization**
\\[ \boldsymbol{U} = \boldsymbol{S}^{-1/2} \\]
This matrix is symmetric (equal to its transpose), so the transformation is simply
\\[\boldsymbol{U}^T \boldsymbol{S} \boldsymbol{U} = \boldsymbol{S}^{-1/2} \boldsymbol{S} \boldsymbol{S}^{-1/2} = \boldsymbol{I} \\]
The advantage of the symmetric orthogonalization is that it is simpler and the new basis is as "close" to the original as possible. We will use the symmetric orthogonalization here.

In order to apply this transformation to the equation \\(\boldsymbol{F}\boldsymbol{C} = \boldsymbol{S}\boldsymbol{C}\boldsymbol{\epsilon}  \\), we need to transform both the Fock matrix and the overlap matrix, since they both depend on the AO basis. That is, we can't just change the basis of one side of the equation, we need to do it to both sides. Mathematically, we need to do this: \\(\boldsymbol{S}^{-1/2}\boldsymbol{F}\boldsymbol{S}^{-1/2} \\) and this: \\(\boldsymbol{S}^{-1/2}\boldsymbol{S} \boldsymbol{S}^{-1/2} \\). This is achieved by multiplying both sides from the left by \\(\boldsymbol{S}^{-1/2}\\) and inserting the identity matrix \\(\boldsymbol{I} = \boldsymbol{S}^{-1/2}\boldsymbol{S}^{1/2}\\) between  \\(\boldsymbol{F}\\) and \\(\boldsymbol{C}\\) on the left and  \\(\boldsymbol{S}\\) and \\(\boldsymbol{C}\\) on the right:

\\[\boldsymbol{F}\boldsymbol{C} = \boldsymbol{S}\boldsymbol{C}\boldsymbol{\epsilon}  \\]

\\[\boldsymbol{S}^{-1/2} \boldsymbol{F} \boldsymbol{S}^{-1/2} \boldsymbol{S}^{1/2} \boldsymbol{C} = \boldsymbol{S}^{-1/2} \boldsymbol{S} \boldsymbol{S}^{-1/2}  \boldsymbol{S}^{1/2}  \boldsymbol{C} \boldsymbol{\epsilon}  \\]

\\[\tilde{\boldsymbol{F}} \tilde{\boldsymbol{C}} = \tilde{\boldsymbol{C}} \boldsymbol{\epsilon} \\]

\\[ \tilde{\boldsymbol{F}} = \boldsymbol{S}^{-1/2}\boldsymbol{F}\boldsymbol{S}^{-1/2} \quad \quad \quad \tilde{\boldsymbol{C}} = \boldsymbol{S}^{1/2}\boldsymbol{C} \\]

We can now easily diagonalize \\(\tilde{\boldsymbol{F}} \\) and transform \\(\tilde{\boldsymbol{C}} \\) back to the original unorthogonal basis \\(\boldsymbol{C} =\boldsymbol{S}^{-1/2}\tilde{\boldsymbol{C}}\\), build a new density matrix, Fock matrix, get a new energy, transform the Fock matrix again, diagonalize... etc. etc. until self-consistency is achieved in the LCAO-MO coefficient matrix \\(\boldsymbol{C}\\) (and thus the density matrix).


The energy expression is found by taking the restricted case of the first Slater-Condon rule

\\[E = 2 \sum\limits_{i}^{N/2} \langle \phi_i | \hat{h}(i) | \phi_i \rangle +  \sum\limits_{i}^{N/2} \sum\limits_{j}^{N/2} 2 \langle \phi_i \phi_j | \hat{g}(i,j) | \phi_i \phi_j \rangle - \langle \phi_i \phi_j | \hat{g}(i,j) | \phi_j \phi_i \rangle\\]

and expanding each spatial orbital in a basis. This is left as an  excercise. The final result is:
\\[E = 2 \sum\limits_{pq} D_{pq} \langle \chi_p \mid \hat{h} \mid \chi_q \rangle  + \sum\limits_{pqrs} D_{pq} D_{rs} [ 2 \langle \chi_p \chi_r \mid \hat{g} \mid \chi_q \chi_s \rangle - \langle \chi_p \chi_r \mid \hat{g} \mid \chi_s \chi_q \rangle ] \\]

Of course, the nuclear repulsion energy would be added to this as a constant.


So, one algorithm for solving the Roothaan-Hall equations and obtaining whats often referred to as the "Restricted Hartree-Fock" (RHF) energy is:

1. Collect all one and two-electron integrals in a matrix, form the orthogonalizer \\(\boldsymbol{S}^{-1/2} \\) 

2. Guess \\(\boldsymbol{D} = \boldsymbol{0}\\)

3. Build \\(\boldsymbol{F}\\)

4. Compute the energy 

5. Diagonalize \\(\tilde{\boldsymbol{F}} = \boldsymbol{S}^{-1/2}\boldsymbol{F}\boldsymbol{S}^{-1/2} \\) to get \\(\tilde{\boldsymbol{C}}\\) and \\(\boldsymbol{\epsilon}\\)

6. Backtransform to unorthogonalized AO basis \\(\boldsymbol{C} =\boldsymbol{S}^{-1/2}\tilde{\boldsymbol{C}}\\)

7. Compute the new density matrix \\(\boldsymbol{D} = \sum\limits_j^{N/2} C_{rj}^* C_{sj}\\)

8. If new \\(\boldsymbol{D}\\) and old \\(\boldsymbol{D}\\) differ by too much, return to step 3.


### Epilogue 1: Where is the Slater determinant?

Taking a look at our RHF algorithm above... Where is our wavefunction? Where is the _Slater determinant_? Didn't we go through all this fuss about representing our wavefunction as a Slater determinant, but it doesn't come up anywhere in solving the Roothaan-Hall equations. It doesn't come up anywhere in our Python code. You're right! It's never actually directly used. The  derivation of the first Slater-Condon rule as well as the Hartree-Fock equations (and therefore, the Roothaan-Hall equations) all assumed a Slater determinant form of the wavefunction. So its never used directly in the algorithm, but it influenced the form of our equations. We _could_ build a Slater determinant wavefunction from our LCAO-MO coefficients combined with our AO basis functions, but it is not needed for computing the energy.

### Epilogue 2: Two-electron integral notations

We presented everything above in terms of _physicist's notation_ (Dirac notation) for two electron integrals:

\\[ \langle \chi_i \chi_j \mid \hat{g} \mid \chi_k \chi_l \rangle = \int \chi_i^*(\boldsymbol{r}_1)\chi_j^*(\boldsymbol{r}_2) \frac{1}{\boldsymbol{r}_{12}}\chi_k(\boldsymbol{r}_1)\chi_l(\boldsymbol{r}_2) d\boldsymbol{r}_1 d\boldsymbol{r}_2 \\]

where we group the real and complex conjugate functions together.

However, in a lot of quantum chemistry literature, some prefer _chemist's notation_ (Pople notation), where the functions are ordered by which electron's coordinates they correspond to. These are often represented with square brackets or parathenses:

\\[ [\chi_i \chi_j \mid \hat{g} \mid \chi_k \chi_l ] = \int \chi_i^*(\boldsymbol{r}_1)\chi_j(\boldsymbol{r}_1) \frac{1}{\boldsymbol{r}_{12}}\chi_k^*(\boldsymbol{r}_2)\chi_l(\boldsymbol{r}_2) d\boldsymbol{r}_1 d\boldsymbol{r}_2 \\]



**Psi4 gives you two-electron integrals in chemist's notation. Therefore you need to transpose Psi4's two-electron integral tensor if you are referencing physicist's notation equations when coding. Alternatively, you can translate all equations into chemists notation and use Psi4's two-electron integrals as given.**

### Epilogue 3: Is our energy computation correct?
You may or may not have noticed that we are computing our energy expression in the _unorthogonalized AO basis_. This may sound in alarm in your head: "Wait a second, we derived the energy expression under the assumption our orbitals are orthonormal, but we are computing the energy with non-orthonormal atomic orbitals!" You'd be right. Here is our energy expression:

\\[E = 2 \sum\limits_{pq} D_{pq} \langle \chi_p \mid \hat{h} \mid \chi_q \rangle  + \sum\limits_{pqrs} D_{pq} D_{rs} [ 2 \langle \chi_p \chi_r \mid \hat{g} \mid \chi_q \chi_s \rangle - \langle \chi_p \chi_r \mid \hat{g} \mid \chi_s \chi_q \rangle ] \\]

All of those \\(\chi \\)'s are AO's, a set of _non-orthonormal_ AO's. This can be seen most easily by looking at your overlap matrix in this AO basis: it is not an identity matrix. So, are we disobeying our assumption that the MO's (not AO's) are orthonormal, which was needed to derive the 1st Slater-Condon rule, hence our energy expression above?

First, we note the following is true for our MO coefficient matrix \\(\boldsymbol{C}\\): it can transform our AO basis to an orthonormal AO basis. You can verify this empirically if you like:

\\[\boldsymbol{C}^T \boldsymbol{S} \boldsymbol{C} = \boldsymbol{1}  \\]

Why is this true? Consider the definition of \\(\boldsymbol{C}\\):

\\[\boldsymbol{C} = \boldsymbol{S}^{-1/2} \boldsymbol{\tilde{C}} \\]

taking the transpose of both sides (\\(\boldsymbol{S}^{-1/2} \\) is symmetric):

\\[\boldsymbol{C}^T = \boldsymbol{\tilde{C}}^T \boldsymbol{S}^{-1/2}  \\]

So we have that 

\\[\boldsymbol{C}^T \boldsymbol{S} \boldsymbol{C} =\boldsymbol{\tilde{C}}^T \boldsymbol{S}^{-1/2} \boldsymbol{S} \boldsymbol{S}^{-1/2} \boldsymbol{\tilde{C}} = \boldsymbol{\tilde{C}}^T  \boldsymbol{\tilde{C}} = \boldsymbol{1} \\]

The last equality is true because the eigenvectors of a symmetric matrix are _always orthonormal_, and \\(\boldsymbol{\tilde{C}}\\) is a matrix of eigenvectors of \\(\boldsymbol{\tilde{F}}\\), which is a symmetric matrix. You can verify this empirically using your RHF code.

Writing our derived expression in summation notation,

\\[ \boldsymbol{C}^T \boldsymbol{S} \boldsymbol{C} = \sum\limits_{pq} C_{pr} S_{pq} C_{qs} = \boldsymbol{1}_{rs}\\]
\\[ \sum\limits_{pq} C_{pr} \langle \chi_p \mid \chi_q \rangle C_{qs} = \boldsymbol{1}_{rs}\\]

The above can be re-expressed using our definition for the AO basis expansion of MO's:

\\[ \langle \phi_r \mid \phi_s \rangle = \boldsymbol{1}_{rs} \\]

What this tells us is the overlaps of all of our molecular orbitals form an identity matrix. Thus, our MO coefficients \\(\boldsymbol{C}\\) (_not_ \\(\boldsymbol{\tilde{C}}\\)) have the effect of transforming our AO-basis overlap matrix to the orthogonalized MO basis when the above operation is performed. So, we might write the following:

\\[\boldsymbol{C}^T \boldsymbol{S}^{AO} \boldsymbol{C} = \boldsymbol{S}^{MO} \\]
\\[\sum\limits_{pq} C_{pr} S^{AO}_{pq} C_{qs} = \boldsymbol{S}^{MO}_{rs} \\]

Summing both sides over the just the \\(n\\) doubly occupied MO's we obtain

\\[\sum\limits_{pq} D_{pq} S_{pq} = n \\]

So, a sum over the occupied dimesions of our MO-basis overlap (an identity matrix) gives the same result as a full contraction of our AO-basis overlap and AO-basis density matrix. It's a bit more involved, but one can show the same is true for our energy contributions in our energy expression: contractions of our AO-basis density matrix with the integral arrays constructed from non-orthogonal AO basis functions gives the same result as using the _orthonormal_ MO-basis quantities directly.
Similar to the above relationships between the AO-basis overlap and MO-basis overlap, the one-electron integrals obey the following:

\\[ \boldsymbol{H}^{MO} = \boldsymbol{C}^T \boldsymbol{H} \boldsymbol{C}   \\]

In more explicit notation the above really means:

\\[\boldsymbol{H}^{MO} = \langle \phi_i \mid \hat{h} \mid \phi_j \rangle = \sum\limits_{pq} C_{pi} \langle \chi_p \mid \hat{h} \mid \chi_q \rangle C_{qj} \\]

We can write a similar transformation for the two electron integrals:

\\[\boldsymbol{G}_{ijkl}^{MO} = \langle \phi_i \phi_j \mid \hat{g} \mid \phi_k \phi_l \rangle =  \sum\limits_{pqrs} C_{pi}C_{qj}C_{rk}C_{sl} \langle \chi_p \chi_q \mid \hat{g} \mid \chi_r \chi_s \rangle \\]

As for the density matrix, 
\\[\boldsymbol{D}^{MO} = \boldsymbol{C}^{-1} \boldsymbol{D}^{AO} (\boldsymbol{C}^{-1})^{T}   \\]

So, the energy expression in terms of our MO basis quantities is:

\\[E = 2 \sum\limits_{pq} D^{MO}_{pq} H^{MO}_{pq} + \sum\limits_{pqrs} D_{pq}^{MO} D_{rs}^{MO} [ 2 G^{MO}_{prqs} - G^{MO}_{prsq} ] \\]

The above energy expression is equivalent to our AO-basis energy expression,

\\[E = 2 \sum\limits_{pq} D_{pq} \langle \chi_p \mid \hat{h} \mid \chi_q \rangle  + \sum\limits_{pqrs} D_{pq} D_{rs} [ 2 \langle \chi_p \chi_r \mid \hat{g} \mid \chi_q \chi_s \rangle - \langle \chi_p \chi_r \mid \hat{g} \mid \chi_s \chi_q \rangle ] \\]

These two energy expressions are not equivalent in a _term-by-term_ sense (if you write out each term in the sums, they are not the same for the AO basis computation and MO basis computation).  However, the contraction of the AO density matrix with the AO-basis integral quantities each evaluate to the same total energy as using integrals over orthogonal MO's. 

There's yet another interesting way to get the same RHF energy. If we **really wanted to,** we could actually use \\(\tilde{\boldsymbol{C}} \\) directly in our energy computation. Recall we obtained \\(\tilde{\boldsymbol{C}} \\) by changing our AO basis with the 'orthogonalizer' \\(\boldsymbol{S}^{-1/2} \\), and we had to transform \\(\boldsymbol{F}\\) and \\(\boldsymbol{S}\\) to \\(\tilde{\boldsymbol{F}} = \boldsymbol{S}^{-1/2}\boldsymbol{F}\boldsymbol{S}^{-1/2} \\) and \\(\boldsymbol{I} = \boldsymbol{S}^{-1/2}\boldsymbol{S}\boldsymbol{S}^{-1/2} \\). To use \\(\tilde{\boldsymbol{C}} \\) directly, we just build our density matrix directly from it:

\\[ \tilde{\boldsymbol{D}}_{pq} = \sum\limits_j^{N/2} \tilde{C}^*_{pj} \tilde{C}_{qj} \\]

We need to also transform our one and two electron integrals with the same \\(\boldsymbol{S}^{-1/2} \\) transformation:
\\[\tilde{\boldsymbol{H}} = \boldsymbol{S}^{-1/2} \boldsymbol{H} \boldsymbol{S}^{-1/2}  \\]
\\[\tilde{\boldsymbol{G}}_{ijkl} = \sum\limits_{pqrs} S^{-1/2}_{pi}S^{-1/2}_{qj}S^{-1/2}_{rk}S^{-1/2}_{sl} G_{pqrs} \\]

Our energy expression is just:

\\[E = 2 \sum\limits_{pq} \tilde{D}_{pq} \tilde{H}_{pq} + \sum\limits_{pqrs} \tilde{D}_{pq} \tilde{D}_{rs} [ 2 \tilde{G}_{prqs} - \tilde{G}_{prsq} ] \\]

Why does this energy expression give the same energy as the AO-basis and MO-basis energy expressions given earlier? We are transforming every integral quantity in terms of non-orthogonal AO basis functions \\(\chi\\) to a basis of _orthogonal_ AO basis functions \\(\tilde{\chi}\\):

\\[\boldsymbol{I} = \boldsymbol{S}^{-1/2} \boldsymbol{S} \boldsymbol{S}^{-1/2}\\]

\\[\boldsymbol{I} = \langle \tilde{\chi}_r \mid \tilde{\chi}_s \rangle = \sum\limits_{pq} S^{-1/2}_{rp} \langle \chi_p \mid \chi_q \rangle S^{-1/2}_{sq}  \\]

So, our AO basis functions are orthonormal, and our coefficients \\(\tilde{\boldsymbol{C}} \\) are also orthonormal, since they are the eigenvectors of a symmetric matrix (\\(\tilde{\boldsymbol{F}} \\)). It follows that molecular orbitals constructed from these orthonormal AO basis functions and orthonormal coefficient eigenvectors are also orthonormal. So the above energy expression is in terms of orthonormal MO's.


The conclusion of all of this is simple: our application of the 1st Slater-Condon rule energy expression **is correct**, because the molecular orbitals **are orthonormal**. Perhaps counter-intuitively, our set of non-orthogonal MO expansion coefficients \\(\boldsymbol{C}\\) and non-orthogonal AO's \\(\chi\\) **construct orthogonal MO's** (concluded from our first investigation showing \\(\boldsymbol{C}^T \boldsymbol{S} \boldsymbol{C} = \boldsymbol{1} \\)). Thus, the energy expression in terms of these quantities is correct. 

We can also compute the energy using MO-basis integral arrays \\(\boldsymbol{H}^{MO}\\) and \\( \boldsymbol{G}^{MO} \\). This MO-basis energy expression is also in terms of orthogonal MO's. 

Finally, we showed how one could use the orthonormal expansion coefficients \\(\tilde{\boldsymbol{C}} \\) combined with orthonormal AO basis functions \\(\tilde{\chi}\\). Together these construct orthonormal MO's, so transforming every integral quantity appropriately by \\(\boldsymbol{S}^{-1/2}\\) in our energy expression gives us the correct result.

Our use of the 1st-Slater-Condon rule energy expression is valid, as long as we are consistent with our AO basis and expansion coefficients across each term in our energy expression, and these AO basis functions and expansion coefficients construct orthogonal MO's. 


