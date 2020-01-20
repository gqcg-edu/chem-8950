## Closed-Shell Case of the First Slater-Condon Rule
Previously we derived that applying the time-independent Schrodinger equation given 1. a Slater determinant wavefunction made up of a linear combination of products of one-electron wavefunctions (spin orbitals), and 2. our electronic Hamiltonian gives the following energy expression

\\[E =  \sum\limits_{i}^N \langle\psi_i^i|\hat{h}(i)|\psi_i^i \rangle + \sum\limits_{i<j}^N  \langle \psi_i^i\psi_j^j|\hat{g}(i,j)|\psi_i^i\psi_j^j \rangle - \langle \psi_i^i\psi_j^j|\hat{g}(i,j)|\psi_j^i\psi_i^j \rangle \\]

Recall that every \\(\psi_p^q\\) is a product of spatial and spin functions \\(\phi_p(\boldsymbol{r}^q)\omega(s^q) \\). For our purposes, the spin function \\(\omega \\) always takes on one of two flavors, which we denote \\( \alpha \\) and \\(\beta\\); spin up and spin down. We assumed previously when deriving the energy expression the following: each spatial orbital is not necessarily the same as any other spatial orbital. However, we _could_ put two electrons, of opposite spin, in the exact same spatial orbital. This is indeed what Nature seems to do. So, we can explore what happens to the energy expression above if we have \\(N = 2n\\) electrons in \\(n\\) doubly-occupied spatial orbitals (where the spin functions \\(\omega\\) of each of electron in the spatial orbital are \\(\alpha\\) and \\(\beta\\)). 

### One-electron term
The summation over \\(N\\) electrons can be converted to a summation over \\(MO\\) doubly-occupied molecular orbitals and a sum over two possible spin functions, \\(\alpha\\) or \\(\beta\\). When we do this, the electrons lose their own coordinate designation; instead they are labeled by which MO and which spin state they are in.

\begin{align}
\sum\limits_{i}^N \langle\psi_i^i|\hat{h}(i)|\psi_i^i \rangle &= \sum\limits_{i}^{MO} \sum\limits_{\omega=\alpha,\beta} \langle\phi_i \omega|\hat{h}(i)|\phi_i \omega \rangle \\
&= \sum\limits_{i}^{MO} \sum\limits_{\omega=\alpha,\beta}  \langle\phi_i |\hat{h}(i)|\phi_i  \rangle \langle \omega | \omega \rangle \\
&= \sum\limits_{i}^{MO}\langle\phi_i |\hat{h}(i)|\phi_i  \rangle \langle \alpha | \alpha \rangle + \langle\phi_i |\hat{h}(i)|\phi_i  \rangle \langle \beta | \beta \rangle
\end{align}

Spin functions are orthonormal, so the self-overlaps are equal to one. Each doubly-occupied spatial orbital contributes _two equivalent amounts of energy_, one for the \\(\alpha\\) spin electron and one for the \\(\beta\\) spin electron. These two units of energy are equivalent, so we can just write

\\[\sum\limits_{i}^N \langle\psi_i^i|\hat{h}(i)|\psi_i^i \rangle = 2\sum\limits_{i}^{MO}\langle\phi_i |\hat{h}(i)|\phi_i  \rangle \\]


### Two-electron term
\\[ \sum\limits_{i<j}^N  \langle \psi_i^i\psi_j^j|\hat{g}(i,j)|\psi_i^i\psi_j^j \rangle - \langle \psi_i^i\psi_j^j|\hat{g}(i,j)|\psi_j^i\psi_i^j \rangle \\]
To make things a little clearer later on, we can re-express the sum over \\(i<j\\) as

\\[\sum\limits_{i<j}^N \rightarrow \frac{1}{2}\sum\limits_{i}^N \sum\limits_{j}^N\\]
To see why,  you can think of an NxN matrix, label one axis \\(i\\) and the other axis \\(j\\). Noting that two-electron integrals cancel each other when \\(i=j\\), and the \\(\frac{1}{2}\\) takes care of double-counting, its apparent that each summation expression covers all entries in the upper (or lower) triangle of this matrix, not including the diagonal. So, we now write this in terms of our new summation:

\begin{align}
\frac{1}{2}\sum\limits_{i}^N \sum\limits_{j}^N \langle \psi_i^i\psi_j^j|\hat{g}(i,j)|\psi_i^i\psi_j^j \rangle - \langle \psi_i^i\psi_j^j|\hat{g}(i,j)|\psi_j^i\psi_i^j \rangle
\end{align}

Let's look at just the 1st (Coulomb) term. Again, we can convert the sums over N electrons into a sum over the number of MO's, but this time we have to be careful to distinguish (label) the spin functions for electrons in MO \\(i\\) and electrons in MO \\(j\\).  
\begin{align}
\frac{1}{2}\sum\limits_{i}^N \sum\limits_{j}^N \langle \psi_i^i\psi_j^j|\hat{g}(i,j)|\psi_i^i\psi_j^j \rangle &= 
\frac{1}{2}\sum\limits_{i}^{MO} \sum\limits_{\omega_i = \alpha,\beta}  \sum\limits_{j}^{MO}  \sum\limits_{\omega_j = \alpha,\beta} \langle \phi_i \omega_i,\phi_j \omega_j | \hat{g}(i,j) | \phi_i \omega_i, \phi_j \omega_j \rangle  \\
&= \frac{1}{2}\sum\limits_{i}^{MO} \sum\limits_{\omega_i = \alpha,\beta}  \sum\limits_{j}^{MO}  \sum\limits_{\omega_j = \alpha,\beta} \langle \phi_i \phi_j | \hat{g}(i,j) | \phi_i \phi_j \rangle \langle \omega_i \omega_j | \omega_i \omega_j \rangle\\
\end{align}


What does the "spin integral" \\(\langle \omega_i \omega_j | \omega_i \omega_j \rangle \\) **do**? Well, first recall our convention for writing all of these two electron integrals: the electron coordinates are kept in order from left to right in the bra and ket \\(\langle 1,2|1,2 \rangle \\). Since there is no operator coupling the two electrons in the case of the spin integral, we can separate it:
\\[ \int \omega_i^*(s_1) \omega_j^*(s_2) \omega_i(s_1) \omega_j(s_2) ds_1 ds_2 = \int \omega_i^*(s_1)  \omega_i(s_1) ds_1 \int \omega_j^*(s_2) \omega_j(s_2) ds_2 \\]

Expanding the spin function summations gives
\begin{align}
\sum\limits_{\omega_i = \alpha,\beta}\sum\limits_{\omega_j = \alpha,\beta} \langle \omega_i \omega_j | \omega_i \omega_j \rangle &=  
\sum\limits_{\omega_i = \alpha,\beta}\sum\limits_{\omega_j = \alpha,\beta} \langle \omega_i  | \omega_i  \rangle \langle \omega_j  | \omega_j  \rangle \\ 
&= \langle \alpha | \alpha \rangle \langle \alpha | \alpha \rangle + \langle \alpha | \alpha \rangle \langle \beta | \beta \rangle + \langle \beta | \beta \rangle \langle \alpha | \alpha \rangle + \langle \beta | \beta \rangle \langle \beta | \beta \rangle \\
&= 4
\end{align}

So, every coulomb two electron integral \\(\langle \phi_i\phi_j  | \hat{g}(i,j) | \phi_i \phi_j  \rangle \\) gets multiplied by a factor of 4 due to the spin integration:

\\[\frac{1}{2}\sum\limits_{i}^{MO} \sum\limits_{\omega_i = \alpha,\beta}  \sum\limits_{j}^{MO}  \sum\limits_{\omega_j = \alpha,\beta} \langle \phi_i \phi_j | \hat{g}(i,j) | \phi_i \phi_j \rangle \langle \omega_i \omega_j | \omega_i \omega_j \rangle = \frac{1}{2}\sum\limits_{i}^{MO}\sum\limits_{j}^{MO}  4 \langle \phi_i \phi_j | \hat{g}(i,j) | \phi_i \phi_j \rangle\\]

For the exchange part, the only difference is that the ket has the two spin orbitals flipped around:
\begin{align}
\frac{1}{2}\sum\limits_{i}^N \sum\limits_{j}^N \langle \psi_i^i\psi_j^j|\hat{g}(i,j)|\psi_j^i\psi_i^j \rangle &= 
\frac{1}{2}\sum\limits_{i}^{MO} \sum\limits_{\omega_i = \alpha,\beta}  \sum\limits_{j}^{MO}  \sum\limits_{\omega_j = \alpha,\beta} \langle \phi_i \omega_i,\phi_j \omega_j | \hat{g}(i,j) | \phi_j \omega_j, \phi_i \omega_i \rangle  \\
&= \frac{1}{2}\sum\limits_{i}^{MO} \sum\limits_{\omega_i = \alpha,\beta}  \sum\limits_{j}^{MO}  \sum\limits_{\omega_j = \alpha,\beta} \langle \phi_i \phi_j | \hat{g}(i,j) | \phi_i \phi_j \rangle \langle \omega_i \omega_j | \omega_j \omega_i \rangle\\
\end{align}

Does this make any difference from before? **Yes**. 
\begin{align}
\sum\limits_{\omega_i = \alpha,\beta}\sum\limits_{\omega_j = \alpha,\beta} \langle \omega_i \omega_j | \omega_j \omega_i \rangle &=  
\sum\limits_{\omega_i = \alpha,\beta}\sum\limits_{\omega_j = \alpha,\beta} \langle \omega_i  | \omega_j  \rangle \langle \omega_j  | \omega_i  \rangle \\ 
&= \langle \alpha | \alpha \rangle \langle \alpha | \alpha \rangle + \langle \alpha | \beta \rangle \langle \beta | \alpha \rangle + \langle \beta | \alpha \rangle \langle \alpha | \beta \rangle + \langle \beta | \beta \rangle \langle \beta | \beta \rangle \\
&= 2
\end{align}

Every exchange integral only gets multiplied by a factor of 2
\\[\frac{1}{2}\sum\limits_{i}^{MO} \sum\limits_{\omega_i = \alpha,\beta}  \sum\limits_{j}^{MO}  \sum\limits_{\omega_j = \alpha,\beta} \langle \phi_i \phi_j | \hat{g}(i,j) | \phi_j \phi_i \rangle \langle \omega_i \omega_j | \omega_j \omega_i \rangle = \frac{1}{2}\sum\limits_{i}^{MO}\sum\limits_{j}^{MO}  2 \langle \phi_i \phi_j | \hat{g}(i,j) | \phi_j \phi_i \rangle\\]

Combining the one electron, coulomb and exchange terms together and cancelling the factor of \\(\frac{1}{2}\\) out front, we have the energy in the special case of a closed-shell Slater determinant:

\\[E = 2 \sum\limits_{i}^{MO} \langle \phi_i | \hat{h}(i) | \phi_i \rangle +  \sum\limits_{i}^{MO} \sum\limits_{j}^{MO} 2 \langle \phi_i \phi_j | \hat{g}(i,j) | \phi_i \phi_j \rangle - \langle \phi_i \phi_j | \hat{g}(i,j) | \phi_j \phi_i \rangle\\]

This is the energy expression used for Restricted Hartree-Fock (RHF). 

