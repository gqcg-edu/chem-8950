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
is a standard eigenvalue problem if \\(\boldsymbol{S} \\) is an identity matrix. If this is the case, the Roothaan-Hall equations are easy to solve; just diagonalize \\(\boldsymbol{F}\\)! Also, if \\(\boldsymbol{S} \\) is an identity matrix, our basis functions are orthonormal. Recall that our energy expression derivation and the Hartree-Fock equations derivation depended on the basis functions being orthonormal. So 'orthogonalizing' our basis is not only is _sufficient_ to make the Roothaan-Hall equations easy to solve, it's also _necessary_ because our equation derivations depended on orbital orthonormality. So, we seek some tranformation \\(\boldsymbol{U}\\) such that 

\\[\boldsymbol{U}^T\boldsymbol{S}\boldsymbol{U} = \boldsymbol{I} \\]

There's two such transformations \\(\boldsymbol{U} \\) that are often used. The first is known as the **canonical orthogonalization**, which is found by diagonalizing the overlap matrix:

\\[\boldsymbol{V}^T\boldsymbol{S}\boldsymbol{V} = \boldsymbol{\Lambda} \\]

\\[\boldsymbol{U} = \boldsymbol{V} \boldsymbol{\Lambda}^{-1/2} \\]

\\[\boldsymbol{U}^T \boldsymbol{S} \boldsymbol{U} = \boldsymbol{\Lambda}^{-1/2} \boldsymbol{V}^T \boldsymbol{S} \boldsymbol{V} \boldsymbol{\Lambda}^{-1/2} = \boldsymbol{\Lambda}^{-1/2} \boldsymbol{\Lambda} \boldsymbol{\Lambda}^{-1/2} = \boldsymbol{I}  \\]
The canonical orthogonalization carries the advantage that we can identify linear dependency in the basis set (an eigenvalue in \\(\boldsymbol{\Lambda}\\) would be near zero!), and **remove** the linear dependency by just setting the corresponding column of \\( \boldsymbol{V}\\) by setting its eigenvalue partner to zero.

The second \\(\boldsymbol{U}\\) used often to orthogonalize our basis is known as the **symmetric orthogonalization**
\\[ \boldsymbol{U} = \boldsymbol{S}^{-1/2} \\]
This matrix is symmetric (equal to its transpose), so the transformation is simply
\\[\boldsymbol{U}^T \boldsymbol{S} \boldsymbol{U} = \boldsymbol{S}^{-1/2} \boldsymbol{S} \boldsymbol{S}^{-1/2} = \boldsymbol{I} \\]
The advantage of the symmetric orthonalization is that it is simpler and the new basis is as "close" to the original as possible. For the rest of this document, we will use the symmetric orthogonalization.

In order to apply this transformation to the equation \\(\boldsymbol{F}\boldsymbol{C} = \boldsymbol{S}\boldsymbol{C}\boldsymbol{\epsilon}  \\), we need to transform both the Fock matrix and the overlap matrix, since they both depend on the AO basis. That is, we can't just change the basis of one side of the equation, we need to do it to both sides. Mathematically, we need to do this: \\(\boldsymbol{S}^{-1/2}\boldsymbol{F}\boldsymbol{S}^{-1/2} \\) and this: \\(\boldsymbol{S}^{-1/2}\boldsymbol{S} \boldsymbol{S}^{-1/2} \\). This is achieved by multiplying both sides from the left by \\(\boldsymbol{S}^{-1/2}\\) and inserting the identity matrix \\(\boldsymbol{I} = \boldsymbol{S}^{-1/2}\boldsymbol{S}^{1/2}\\) between  \\(\boldsymbol{F}\\) and \\(\boldsymbol{C}\\) on the left and  \\(\boldsymbol{S}\\) and \\(\boldsymbol{C}\\) on the right:

\\[\boldsymbol{F}\boldsymbol{C} = \boldsymbol{S}\boldsymbol{C}\boldsymbol{\epsilon}  \\]

\\[\boldsymbol{S}^{-1/2} \boldsymbol{F} \boldsymbol{S}^{-1/2} \boldsymbol{S}^{1/2} \boldsymbol{C} = \boldsymbol{S}^{-1/2} \boldsymbol{S} \boldsymbol{S}^{-1/2}  \boldsymbol{S}^{1/2}  \boldsymbol{C} \boldsymbol{\epsilon}  \\]

\\[ \tilde{\boldsymbol{F}} \tilde{\boldsymbol{C}} = \tilde{\boldsymbol{C}} \boldsymbol{\epsilon} \quad \quad \quad \tilde{\boldsymbol{F}} = \boldsymbol{S}^{-1/2}\boldsymbol{F}\boldsymbol{S}^{1/2} \quad \quad \quad \tilde{\boldsymbol{C}} = \boldsymbol{S}^{1/2}\boldsymbol{C} \\]

We can now easily diagonalize \\(\tilde{\boldsymbol{F}} \\) and transform \\(\tilde{\boldsymbol{C}} \\) back to the original unorthogonal basis \\(\boldsymbol{C} =\boldsymbol{S}^{-1/2}\tilde{\boldsymbol{C}}\\), build a new density matrix, Fock matrix, get a new energy, transform the Fock matrix again, diagonalize... etc. etc. until self-consistency is achieved in the LCAO-MO coefficient matrix \\(\boldsymbol{C}\\) (and thus the density matrix).


The energy expression is found by taking the restricted case of the first Slater-Condon rule

\\[E = 2 \sum\limits_{i}^{N/2} \langle \phi_i | \hat{h}(i) | \phi_i \rangle +  \sum\limits_{i}^{N/2} \sum\limits_{j}^{N/2} 2 \langle \phi_i \phi_j | \hat{g}(i,j) | \phi_i \phi_j \rangle - \langle \phi_i \phi_j | \hat{g}(i,j) | \phi_j \phi_i \rangle\\]

and expanding each spatial orbital in a basis. This is left as an  excercise. The final result is:

TODO TODO TODO!!!

So, one algorithm for solving the Roothaan-Hall equations and obtaining the whats often referred to as the "Restricted Hartree-Fock" (RHF) energy is:

1. Collect all one and two-electron integrals in a matrix, form the orthgonalizer \\(\boldsymbol{S}^{1/2} \\) 

2. Guess \\(\boldsymbol{D} = \boldsymbol{0}\\)

3. Build \\(\boldsymbol{F}\\)

4. Compute the energy 

5. Diagonalize \\(\tilde{\boldsymbol{F}} = \boldsymbol{S}^{-1/2}\boldsymbol{F}\boldsymbol{S}^{1/2} \\) to get \\(\tilde{\boldsymbol{C}}\\) and \\(\boldsymbol{\epsilon}\\)

6. Backtransform to unorthogonalized AO basis \\(\boldsymbol{C} =\boldsymbol{S}^{-1/2}\tilde{\boldsymbol{C}}\\)

7. Compute the new density matrix \\(\boldsymbol{D} = \sum\limits_j^{N/2} C_{rj}^* C_{sj}\\)

8. If new \\(\boldsymbol{D}\\) and old \\(\boldsymbol{D}\\) differ by too much, return to step 3.


### Epilogue 1: Is our energy computation correct?

You may or may not have noticed that we are computing our energy expression in the _unorthogonalized AO basis_. This may sound in alarm in your head: "Wait a second, we derived the energy expression under the assumption our orbitals are orthonormal, but we are computing the energy in an un-orthonormal basis!" You'd be right. But, recall from the excercises of a previous set of notes, you proved that if we subject our orbitals to some linear transformation
\\[ \boldsymbol{\phi'} = \boldsymbol{\phi} \boldsymbol{A} \\] we only scale our wavefunction (Slater determinant) by a scalar quantity, 

\\[ \Phi = \Phi \det(\boldsymbol{A}) \\]


the remnant of orthonormality that once was make everything okay. One _could_ just transform everything to the orthogonalized AO basis (one electron integrals, two electron integrals, etc) but this would be needlessly expensive.


### Epilogue 2: Where is the Slater determinant?

Taking a look at our RHF algorithm above... Where is our wavefunction? Where is the _Slater determinant_? Didn't we go through all this fuss about representing our wavefunction as a Slater determinant, but it doesn't come up anywhere in solving the Roothaan-Hall equations. It doesn't come up anywhere in our Python code. You're right! It's never actually directly used. 
(to be continued)


### Epilogue 3: Two-electron integral notations

We presented everything above in terms of _physicist's notation_ (Dirac) for two electron integrals:

However, in a lot of quantum chemistry literature, some prefer _chemist's notation_ (Pople)




