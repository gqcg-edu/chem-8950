# Molecular Integrals and Basis Sets

Many equations which appeared in our discussion of Hartree-Fock theory included various _integrals_ involving atomic orbital basis functions and operators from our electronic Hamiltonian. These are integrals over all space, and they evaluate to a single value; a floating point number. We always refer to them as "integrals", but more precisely, they are the numerical result of the definite integral over all space of the AO functions and operator in the integrand. A few examples of Hartree-Fock-relevant equations involving such integrals:

* The Fock matrix elements are composed of a sum of integrals over AO basis functions (\\(\chi\\)): electron-kinetic energy integrals (with operator  \\(  \hat{T} = \frac{1}{2} \nabla^2 \\)), electron-nuclear attraction integrals (with operator \\(\hat{V} = \sum\limits_A \frac{Z_A}{\mid\boldsymbol{r} - \boldsymbol{R}_A \mid}  \\)), and two electron integrals (with operator \\(\hat{g} = \frac{1}{|\boldsymbol{r}_1 - \boldsymbol{r}_2|} \\))

\\[\langle \chi_p \mid \hat{f}\mid \chi_q \rangle= \langle \chi_p \mid \hat{T} \mid \chi_q \rangle + \langle \chi_p \mid \hat{V} \mid \chi_q \rangle + \sum\limits_r^m \sum\limits_s^m D_{rs} \left[ 2 \langle \chi_p \chi_r \mid \hat{g} \mid \chi_q \chi_s \rangle - \langle \chi_p \chi_r \mid \hat{g} \mid \chi_s \chi_q \rangle \right]\\] 

* The energy expression (with \\(\hat{h} = \hat{T} + \hat{V} \\))

\\[E = 2 \sum\limits_{pq} D_{pq} \langle \chi_p \mid \hat{h} \mid \chi_q \rangle  + \sum\limits_{pqrs} D_{pq} D_{rs} [ 2 \langle \chi_p \chi_r \mid \hat{g} \mid \chi_q \chi_s \rangle - \langle \chi_p \chi_r \mid \hat{g} \mid \chi_s \chi_q \rangle ] \\]

* The overlap matrix \\(\boldsymbol{S}\\), which is used to construct our 'orthogonalizer' \\(\boldsymbol{S}^{-1/2}\\)

\\[S_{pq} = \langle \chi_p \mid \chi_q \rangle \\]


**We have not yet specified what these AO basis functions \\(\chi\\) are.** The choice of AO basis functions ultimately will determine the quality of the electronic wavefunction (which is almost always built from antisymmetric products of LCAO-MO's, which are linear combinations of the AO basis functions). Thus, the quality of the AO basis functions will determine the quality of everything we compute. How do we judge the appropriateness of our choice of basis functions? Well, they will ideally satisfy many qualities:

* First of all, they better be functions that do not give divergent integrals; otherwise our energy would be infinite. 

* Second, the required molecular integrals above should be evaluatable with these basis functions. 

* Third, they better be similar-looking to known one-electron wavefunctions (the Hydrogen atom wavefunctions). Any other choice would be harder to justify. 

* Fourth, our basis functions should be systematically improvable, and we should be able to approach completeness with respect to the space of all one-electron square-integrable functions.

There are several possible choices which satisfy these criteria. We will focus only on the most common choice of basis functions: Gaussian basis functions.

## Gaussian Basis Functions

An unnormalized "primitive" Cartesian Gaussian function \\(\chi(\boldsymbol{r})\\), which is a function of the three spatial coordinates of the electron  \\( \boldsymbol{r} = (x,y,z)  \\) centered at \\(\boldsymbol{A} = (A_x, A_y, A_z) \\) with exponent \\(\alpha\\) is
\\[ \chi_k(\boldsymbol{r}) = (x - A_x)^{a_x} (y - A_y)^{a_y} (z - A_z)^{a_z} \exp[- \alpha_k (\boldsymbol{r} - \boldsymbol{A})^2 ]  \\]

The set of integers \\(\boldsymbol{a} = (a_x, a_y, a_z) \\) are often called the _angular momenta_ of the Gaussian, or, the _quantum numbers_ of the Gaussian. 
The reason for this is that the prefactor \\( (x - A_x)^{a_x} (y - A_y)^{a_y} (z - A_z)^{a_z} \\) gives an angular dependence analogous to the spherical harmonic part of the Hydrogen atom wavefunctions, and the exponential mimics the radial part of the Hydrogen atom wavefunctions. The coordinates of the center \\(\boldsymbol{A} \\) are almost always the Cartesian coordinates of an atom in the molecule. So, the basis functions are _centered_ on the atoms.

The sum of these three integers \\(\boldsymbol{a} = (a_x, a_y, a_z) \\)  is the total angular momentum of the Gaussian. We _label_ the Gaussians according to their angular momentum vector.  An \\(s\\) Gaussian function has \\((a_x, a_y, a_z) = (0,0,0) \\), that is, a total angular momentum of 0. Guassians with total angular momentum 1 are denoted as \\(p\\) functions, 2 are \\(d\\) functions, 3 \\(f\\), 4 \\(g\\), 5 \\(h\\), 6 \\(i\\)...

We also distinguish between different Gaussians with the same total angular momentum. A few examples:

  * \\(p_x \implies \boldsymbol{a} = (1, 0, 0) \\)
  * \\(p_y \implies \boldsymbol{a} = (0, 1, 0) \\)
  * \\(p_z \implies \boldsymbol{a} = (0, 0, 1) \\)
  * \\(d_{xy} \implies \boldsymbol{a} = (1, 1, 0) \\)    
  * \\(f_{xyz} \implies \boldsymbol{a} = (1, 1, 1) \\)
  * \\(g_{xxxx} \implies \boldsymbol{a} = (4, 0, 0) \\)

### Contracted Gaussians

Originally, _Slater functions_ were used as AO basis functions because they get the "cusp conditions" correct; The behavior of Slater functions at the origin properly mimics the Hydrogen-like radial wavefunctions. Gaussian functions, on the other had, have a smooth hump at the origin. 

![basis.gif](attachment:basis.gif)

However, two-electron integrals over Slater functions are basically impossible to solve analytically, while ERIs over Gaussian functions are relatively simple. It was quickly realized that one could fix the cusp conditions of Gaussian basis functions by just taking a linear combination of them until they are 'pointy enough.' These are referred to as _contracted Gaussians_, which are a sum of primitive Gaussian with same angular momentum and center, but different exponents \\(\alpha_k\\):
\\[\chi^c = \sum\limits_k^K D_{k} \chi_{k} \\]

The length of the linear combination is the _degree of contraction_.

## Reading a basis set file
Below is a Psi4 basis set file for STO-3G Oxygen. The numbers on the left are the exponents \\(\alpha \\) and the numbers on the right are contraction coefficients. 

![psi4basis.png](attachment:psi4basis.png)

There are 5 basis functions described here, all of which are contracted. Rounding to the third decimal place, they are:
\begin{align}
\chi_{s_1} &= 0.154(x - A_x)^0(y - A_y)^0(z-A_z)^0 \exp(-0.131(\boldsymbol{r} - \boldsymbol{A})^2 ) \\
           &+ 0.535(x - A_x)^0(y - A_y)^0(z-A_z)^0 \exp(-0.238(\boldsymbol{r} - \boldsymbol{A})^2 ) \\
           &+ 0.445(x - A_x)^0(y - A_y)^0(z-A_z)^0 \exp(-0.644(\boldsymbol{r} - \boldsymbol{A})^2 ) \\
\chi_{s_2} &= -1.0(x - A_x)^0(y - A_y)^0(z-A_z)^0 \exp(-0.503(\boldsymbol{r} - \boldsymbol{A})^2 ) \\
           &+ 0.399(x - A_x)^0(y - A_y)^0(z-A_z)^0 \exp(-0.117(\boldsymbol{r} - \boldsymbol{A})^2 ) \\
           &+ 0.700(x - A_x)^0(y - A_y)^0(z-A_z)^0 \exp(-0.380(\boldsymbol{r} - \boldsymbol{A})^2 ) \\
\chi_{p_x} &= 0.156(x - A_x)^1(y - A_y)^0(z-A_z)^0 \exp(-0.503(\boldsymbol{r} - \boldsymbol{A})^2 ) \\
           &+ 0.608(x - A_x)^1(y - A_y)^0(z-A_z)^0 \exp(-0.117(\boldsymbol{r} - \boldsymbol{A})^2 ) \\
           &+ 0.392(x - A_x)^1(y - A_y)^0(z-A_z)^0 \exp(-0.380(\boldsymbol{r} - \boldsymbol{A})^2 ) \\
\chi_{p_y} &= 0.156(x - A_x)^0(y - A_y)^1(z-A_z)^0 \exp(-0.503(\boldsymbol{r} - \boldsymbol{A})^2 ) \\
           &+ 0.608(x - A_x)^0(y - A_y)^1(z-A_z)^0 \exp(-0.117(\boldsymbol{r} - \boldsymbol{A})^2 ) \\
           &+ 0.392(x - A_x)^0(y - A_y)^1(z-A_z)^0 \exp(-0.380(\boldsymbol{r} - \boldsymbol{A})^2 ) \\
\chi_{p_z} &= 0.156(x - A_x)^0(y - A_y)^0(z-A_z)^1 \exp(-0.503(\boldsymbol{r} - \boldsymbol{A})^2 ) \\
           &+ 0.608(x - A_x)^0(y - A_y)^0(z-A_z)^1 \exp(-0.117(\boldsymbol{r} - \boldsymbol{A})^2 ) \\
           &+ 0.392(x - A_x)^0(y - A_y)^0(z-A_z)^1 \exp(-0.380(\boldsymbol{r} - \boldsymbol{A})^2 ) \\
\end{align}

For maximum clarity, I have written out the angular momentum prefactors even though most are just equal to 1. So, one might say the following about the above basis set for the oxygen atom: there are 5 _basis functions_, each constructed from a _contraction_ of 3 _primitives_ (the contractions have a _degree_ of 3), there are 2 _s_-type Gaussian basis functions with 0 _angular momentum_ and one _p_ function (we just say there's one p function, everyone knows its really 3 functions in 3 different cartesian directions).

What are the centers in the above expression? Well, they would be the nuclear cartesian coordinates of your oxygen atoms. If you have diatomic oxygen, for example, you would have 10 basis functions; the two sets of 5 on each oxygen would be identical except for the center coordinates.

Sometimes basis functions are 'hybridized' and labeled with "SP". This means nothing more than to take the same and orbital exponents and use them to construct \\(s, p_x, p_y, p_z \\) functions. For clarity, I'm using one of the advanced options on basissetexchange.org which allows "uncontract SPDF" which means remove this hybridization notation and just write the basis functions all separately. Note in the above example \\(\chi_{s_2}\\) and the \\(\chi_p\\) functions have the same orbital exponents. Thus, on basissetexhange.org, if you load this basis set with default settings it will lump the second s function and p function together under the label "SP", give just _one_ column of orbital exponents, and two columns of distinct contraction coefficients for the s and p parts. 


## Integrals over Gaussian Basis Functions

We have 4 kinds of integrals to evaluate for Hartree-Fock: overlap, kinetic, potential, and electron repulsion integrals. We need analytic formulas to evaluate each of these integrals when the AO basis functions are Gaussians. S.F. Boys derived analytic expressions for each of these integrals over only \\(s\\)-type Gaussian functions in his seminal 1950 paper, "Electronic Wavefucntions I. A general method of calculation for the stationary states of any molecular system."

### Notation

We will use notation common in the molecular integrals literature. A single Gaussian primitive is typically just represented by a single bold letter, with each quantity within the Gaussian function simply implied to be labeled with the same letter:
\\[\boldsymbol{a} := \chi_a(\boldsymbol{r}) = (x - A_x)^{a_x} (y - A_y)^{a_y} (z - A_z)^{a_z} \exp[- \alpha_a (\boldsymbol{r} - \boldsymbol{A})^2 ]  \\]

\\[\boldsymbol{b} := \chi_b(\boldsymbol{r}) = (x - B_x)^{b_x} (y - B_y)^{b_y} (z - B_z)^{b_z} \exp[- \alpha_b (\boldsymbol{r} - \boldsymbol{B})^2 ]  \\]

We can write a one-electron integral as:

\\[ (\boldsymbol{a} \mid \hat{O} \mid \boldsymbol{b}) = \int \chi_a^*(\boldsymbol{r}_1) \hat{O} \chi_b(\boldsymbol{r}_1) d\boldsymbol{r}_1 \\]

and a two-electron integral as:
\\[ (\boldsymbol{a} \boldsymbol{b}  \mid \boldsymbol{c} \boldsymbol{d}) = \int \chi_a^*(\boldsymbol{r}_1) \chi_b(\boldsymbol{r}_1) r_{12}^{-1}    \chi_c^*(\boldsymbol{r}_2) \chi_d(\boldsymbol{r}_2) d\boldsymbol{r}_1 d\boldsymbol{r}_2\\]

Here it is just implied that these two Gaussians have distinctly labeled centers, angular momentum vectors, coordinates, and exponents.
So, there are a lot of details being hidden with this notation. Also note we are using chemist's notation for the two electron integrals, simply because that is what everyone does in the molecular integrals literature.


### Integrals over \\(s\\) functions

Integrals over just \\(s\\)-type Gaussian functions have clean, analytic forms. Higher-angular momentum integrals are found by very complicated expressions. Fortunately, it has been discovered that all higher angular momentum integrals can be computed exactly from a linear combination of integrals over \\(s\\) functions. Almost every integrals code used today essentially just computes a bunch of \\(s\\) function integrals and combines them together using various _recursion relations_ which effectively promote the angular momentum on each Gaussian integral to the desired value (e.g., \\((s \mid s) \rightarrow (d_{xy} \mid f_{yzz}) \\)). More on this later!

### Overlap integrals over \\(s\\) functions
A normalized \\(s\\) overlap integral over primitive Gaussians is given by

\\[ (\boldsymbol{a}\mid\boldsymbol{b}) = (s \mid s) = N_a N_b \frac{\pi}{\alpha_a + \alpha_b} \exp \left( \frac{-\alpha_a \alpha_b (\boldsymbol{A} - \boldsymbol{B})^2}{\alpha_a + \alpha_b} \right)\\]

The normalization constant for each Gaussian can be found by

\\[N_a = \frac{\left(\frac{2 \alpha_a}{\pi}\right)^{3/4} (4 \alpha_a)^{(a_x + a_y + a_z) /2}}{\sqrt{(2a_x - 1)!!(2a_y - 1)!!(2a_z - 1)!!}}    = \left(\frac{2 \alpha_a}{\pi}\right)^{3/4}\\]

### Kinetic integrals over \\(s\\) functions

The kinetic energy integral can be greatly simplified if expressed in terms of the overlap integral and defining \\( \omega = \frac{\alpha_a \alpha_b}{ \alpha_a + \alpha_b}\\):

\\[ (\boldsymbol{a}\mid \hat{T} \mid \boldsymbol{b}) = (s\mid\hat{T}\mid s) = (s\mid s) \cdot (3 \omega + 2 \omega^2 \cdot -(\boldsymbol{A} - \boldsymbol{B})^2 )\\]

### Potential integrals over \\(s\\) functions

Now things get more complicated. An electron-nuclear potential energy integral requires one to consider the contribution of every single atom's nucleus in the system. Thus, every potential integral depends on the entire set of Cartesian coordinates of the molecule (\\(\boldsymbol{G}_i\\)) and the nuclear charges (\\(Z_i\\)) of each atom \\(i\\). 

First we introduce

\\[\boldsymbol{P} = \frac{\alpha_a \boldsymbol{A} + \alpha_b \boldsymbol{B}}{\alpha_a + \alpha_b} \\]

\\[x_i = (\alpha_a + \alpha_b) (\boldsymbol{P} - \boldsymbol{G}_i)^2    \\]


The potential integral is then
\\[ (\boldsymbol{a}\mid \hat{V} \mid \boldsymbol{b}) = (s\mid \hat{V} \mid s) = N_a N_b \frac{2 \pi}{\alpha_a + \alpha_b} \exp[(\boldsymbol{A} - \boldsymbol{B})^2 \cdot -\omega] \sum\limits_i^{n}  -Z_i F_0(x_i)   \\]

where the sum is over all \\(n\\) nuclei.

The function \\(F_0(x)\\) is known as the \\(0^{th}\\) _Boys function_. It is a family of functions (\\(F_0, F_1, F_2, ...\\)):
\\[ F_\nu(x) =  \int\limits_0^1 t^{2\nu} \exp(-xt^2) dt = \frac{1}{2x^{\nu + \frac{1}{2}}} \cdot \gamma(\nu + \frac{1}{2}, x) \cdot \Gamma(\nu + \frac{1}{2}) \\]

where \\(\gamma\\) is the lower incomplete Gamma function and \\(\Gamma\\) is the Gamma function.

In the special case of \\(F_0(x)\\) we can ignore the above definition entirely and use the very simplified form:

\\[ F_0(x) = \mathrm{erf}(\sqrt x) \frac{\sqrt{\pi}} {2 \sqrt {x}} \\]

The above expression blows up when \\(x\\) is small, so typically for \\(x < 10^{-8} \\) we use the Taylor expansion:

\\[ 1 - \frac{x}{3} + \frac{x^2}{10} - \frac{x^3}{42} \\]

Well, that's as bad as it gets (two electron integrals are easier; no sum over every interaction with every nucleus). It's not _that_ bad though, right?

### Electron repulsion integrals over \\(s\\) functions
For expressing a two-electron integral over \\(s\\) functions it is convenient to define:

\\[\boldsymbol{P} = \frac{\alpha_a \boldsymbol{A} + \alpha_b \boldsymbol{B}}{\alpha_a + \alpha_b} \\]

\\[\boldsymbol{Q} = \frac{\alpha_c \boldsymbol{C} + \alpha_d \boldsymbol{D}}{\alpha_c + \alpha_d} \\]

\\[T = \frac{(\alpha_a + \alpha_b)(\alpha_c + \alpha_d)}{\alpha_a + \alpha_b + \alpha_c + \alpha_d} (\boldsymbol{P} - \boldsymbol{Q})^2 \\]

\\[K_{AB} = 2^{1/2} \frac{\pi^{5/4}}{\alpha_a + \alpha_b} \exp \left[ -\frac{\alpha_a \alpha_b}{ \alpha_a + \alpha_b} (\boldsymbol{A} - \boldsymbol{B})^2 \right] \\]

\\[K_{CD} = 2^{1/2} \frac{\pi^{5/4}}{\alpha_c + \alpha_d} \exp \left[ -\frac{\alpha_c \alpha_d}{ \alpha_c + \alpha_d} (\boldsymbol{C} - \boldsymbol{D})^2 \right] \\]

Our two electron integrals over \\(s\\) functions can now be expressed simply as

\\[(\boldsymbol{ab} \mid \boldsymbol{cd}) = (s s \mid s s) = N_a N_b N_c N_d (\alpha_a + \alpha_b + \alpha_c + \alpha_d )^{-1/2} K_{AB} K_{CD} F_0(T) \\]

where \\(F_0\\) is the same Boys function as defined previously.

### Recap

Some of the above expressions are messy, but at their core they are very simple: every integral is some function of the Gaussian exponents and Gaussian centers, combined with some constants, exponentials, dot products, etc. The only weird part is the Boys function, that only appears in the case of potential integrals and two-electron integrals.


### Higher angular momentum integrals

Primitive Gaussians have the following property when differentiated w.r.t. some Cartesian component of their center:

\\[\frac{\partial}{\partial A_i}(\boldsymbol{a}) = 2 \alpha_a (\boldsymbol{a} + \boldsymbol{1}_i) - a_i (\boldsymbol{a} - \boldsymbol{1}_i) \\]

This formula implies that for all one-electron integrals, the following is true:

\\[\frac{\partial}{\partial A_i} (\boldsymbol{a} \mid \boldsymbol{b}) = 2 \alpha_a ((\boldsymbol{a} + \boldsymbol{1}_i) \mid \boldsymbol{b}) - a_i ((\boldsymbol{a} - \boldsymbol{1}_i) \mid \boldsymbol{b} ) \\]

The expression \\((\boldsymbol{a} + \boldsymbol{1}_i)\\) means "take the Gaussian \\(\boldsymbol{a}\\) and promote the \\(i^{th}\\) cartesian component angular momentum by 1." Mathematically, in the \\( (i - A_i)^{a_i} \\) prefactor of the Gaussian, simply add 1 to \\(a_i\\). In the equation above, the factor \\(a_i ((\boldsymbol{a} - \boldsymbol{1}_i)\mid \boldsymbol{b} ) \\) means to take the **current** angluar momentum on Gaussian \\(\boldsymbol{a}\\) and multiply by the integral \\(((\boldsymbol{a} - \boldsymbol{1}_i) \mid \boldsymbol{b} ) \\), **if that integral exists**. As you might expect, if \\((\boldsymbol{a} \mid \boldsymbol{b} )\\) was an integral over two _s_ functions, there would be no lower angular momentum integral \\(((\boldsymbol{a} - \boldsymbol{1}_i) \mid \boldsymbol{b} ) \\), so this term is zero in that case. There is an analgous expression for two electron integrals:

\\[\frac{\partial}{\partial A_i} (\boldsymbol{a}\boldsymbol{b}\mid \boldsymbol{c}\boldsymbol{d}) = 2 \alpha_a ((\boldsymbol{a} + \boldsymbol{1}_i) \boldsymbol{b} \mid \boldsymbol{c} \boldsymbol{d}) - a_i ((\boldsymbol{a} - \boldsymbol{1}_i) \boldsymbol{b} \mid \boldsymbol{c} \boldsymbol{d}) \\]

Not only do these relations allow us to differentiate our integrals w.r.t. nuclear coordinates (needed for analytic gradient methods!), but they also give an expression for how to get higher angular momentum integrals from lower angular momentum integrals by the simple rearrangement:

\\[((\boldsymbol{a} + \boldsymbol{1}_i) \mid \boldsymbol{b} ) = \frac{1}{2 \alpha_a} \left( \frac{\partial}{\partial A_i} (\boldsymbol{a} \mid \boldsymbol{b}) + a_i ((\boldsymbol{a} - \boldsymbol{1}_i) \mid \boldsymbol{b})  \right)\\]

\\[((\boldsymbol{a} + \boldsymbol{1}_i) \boldsymbol{b} \mid \boldsymbol{c}\boldsymbol{d} ) = \frac{1}{2 \alpha_a} \left( \frac{\partial}{\partial A_i} (\boldsymbol{a}\boldsymbol{b} \mid \boldsymbol{c}\boldsymbol{d}) + a_i ((\boldsymbol{a} - \boldsymbol{1}_i ) \boldsymbol{b}\mid \boldsymbol{c}\boldsymbol{d})  \right)\\]


These derivative relations have been used to derive various **recursion relations**, which describe how to get an integral over arbitrary angular momentum Gaussians in terms of linear combinations of lower angular momentum Gaussians.

### Recurrence Relations

The most famous of all recurrence relations are those of Obara and Saika (OS) given in their paper "Efficient recursive computation of
molecular integrals over Cartesian
Gaussian functions", 1985. They are not the only set of useful recursion relations (see the review by Peter Gill, "Molecular integrals over Gaussian basis functions", 1994). 

The OS recursion was derived from the above derivative relations. We will not derive them here, but instead just present the result and how the recrusions work. We only present the recursion relations (RR) for overlap and two-electron integrals, as the overlap RR is somewhat analogous to the kinetic integral RR, and the two-electron integral RR is analogous to the potential integral RR.

Each of these relations apply to **unnormalized** Gaussians. After the recursion is performed, the appropriate normalization constants are multiplied in, as well as the contraction coefficients, and then the integrals are summed according to the contractions defined by the basis set.

### Overlap Integral Recurrence Relation

\\[ ((\boldsymbol{a} + \boldsymbol{1}_i) \mid \boldsymbol{b}) = (P_i - A_i)(\boldsymbol{a} \mid \boldsymbol{b}) + \frac{1}{2 \alpha_a} a_i ((\boldsymbol{a} - \boldsymbol{1}_i) \mid \boldsymbol{b}) + \frac{1}{2\alpha_a} b_i ( \boldsymbol{a} \mid (\boldsymbol{b} - \boldsymbol{1}_i))\\]

where \\(\boldsymbol{P}\\) is the same as before:
\\[\boldsymbol{P} = \frac{\alpha_a \boldsymbol{A} + \alpha_b \boldsymbol{B}}{\alpha_a + \alpha_b} \\]
and \\(b_i\\) is analogous to the definition of \\(a_i\\): it is the current angular momentum on center \\(\boldsymbol{b}\\) for integral \\((\boldsymbol{a} \mid \boldsymbol{b}) \\).


### Electron Repulsion Integral Recurrence Relation

\\[ ((\boldsymbol{a} + \boldsymbol{1}_i) \boldsymbol{b} \mid \boldsymbol{c} \boldsymbol{d})^{(m)} = (P_i - A_i) (\boldsymbol{a}\boldsymbol{b} \mid \boldsymbol{c}\boldsymbol{d})^{(m)} + (W_i - P_i) (\boldsymbol{a}\boldsymbol{b} \mid \boldsymbol{c}\boldsymbol{d})^{(m+1)}\\]

\\[+ \frac{a_i}{\alpha_a} \left[ ((\boldsymbol{a} - \boldsymbol{1}_i)\boldsymbol{b} \mid \boldsymbol{c}\boldsymbol{d})^{(m)} - \frac{\eta}{\eta + \zeta} ((\boldsymbol{a} - \boldsymbol{1}_i)\boldsymbol{b} \mid \boldsymbol{c}\boldsymbol{d})^{(m+1)}  \right]\\]

\\[+ \frac{b_i}{\alpha_a} \left[ (\boldsymbol{a}(\boldsymbol{b} - \boldsymbol{1}_i) \mid \boldsymbol{c}\boldsymbol{d})^{(m)} - \frac{\eta}{\eta + \zeta} (\boldsymbol{a} (\boldsymbol{b} - \boldsymbol{1}_i)\mid \boldsymbol{c}\boldsymbol{d})^{(m+1)}  \right]\\]

\\[+ \frac{c_i}{2(\zeta + \eta)} ( \boldsymbol{ab}\mid (\boldsymbol{c} - \boldsymbol{1}_i) \boldsymbol{d})^{(m+1)} + \frac{d_i}{2(\zeta + \eta)}( \boldsymbol{ab}\mid  \boldsymbol{c}(\boldsymbol{d} - \boldsymbol{1}_i))^{(m+1)} \\]

Where

\\[\boldsymbol{P} = \frac{\alpha_a \boldsymbol{A} + \alpha_b \boldsymbol{B}}{\alpha_a + \alpha_b} \\]

\\[\boldsymbol{Q} = \frac{\alpha_c \boldsymbol{C} + \alpha_d \boldsymbol{D}}{\alpha_c + \alpha_d} \\]

\\[\zeta = \alpha_a + \alpha_b \\]

\\[\eta = \alpha_c + \alpha_d \\]

\\[\boldsymbol{W} = \frac{\zeta P_i + \eta Q_i}{\zeta + \eta}  \\]


The "m" superscripts denote which Boys function \\(F_m(T) \\) to use in computing the integral.  The "true" integral corresponds to \\(m=0\\).
The terms of each equation can be swapped to obtain expressions for promoting every other center as well. Obara and Saika in their paper list out the analytic forms for all two electron integrals involving _s_ and _p_ functions obtained from this recursion.

![os.png](attachment:os.png)