# Uno

**Contents**
- [Installation](#installation)
- [How to cite Uno](#how-to-cite-uno)
- [Credits](#credits)

Uno (Unifying Nonlinear Optimization) is a C++ framework for solving nonlinearly constrained optimization problems of the form:

$$
\begin{align}
\min_{x \in \mathbb{R}^n}  & ~f(x) \\
\text{s.t.}                & ~c_L \le c(x) \le c_U \\
                           & ~x_L \le x \le x_U \\
\end{align}
$$

where $f: \mathbb{R}^n \rightarrow \mathbb{R}$ and $c: \mathbb{R}^n \rightarrow \mathbb{R}^m$ are (ideally twice) continuously differentiable.

Uno unifies Lagrange-Newton (essentially **SQP** and **interior-point**) methods that iteratively solve the optimality (KKT) conditions with Newton's method. It breaks them down into a set of building blocks that interact with one another. You can combine these strategies in a ton of different ways via options. Uno also implements **presets**, that is strategy combinations that mimic existing solvers:
* `filtersqp` mimics filterSQP (trust-region feasibility restoration filter SQP method with exact Hessian);
* `ipopt` mimics IPOPT (line-search feasibility restoration filter barrier method with exact Hessian and primal-dual inertia correction).

Uno can be used via its [AMPL/nl](interfaces/AMPL/README.md), [Julia](interfaces/Julia/README.md), [Python](interfaces/Python/README.md), [C](interfaces/C/README.md), and [Fortran](interfaces/Fortran/README.md) interfaces.

## Installation

See the [installation guide](INSTALL.md) file for instructions on how to compile Uno from source or use the precompiled libraries and executables.

## Getting started

To get started with Uno, check out the [official documentation](https://unosolver.readthedocs.io).

## How to cite Uno

Our Uno paper was accepted in the Mathematical Programming Computation journal on Feb 22, 2026.
Our preprint is available on [ResearchGate](https://www.researchgate.net/publication/397446552_Implementing_a_unified_solver_for_nonlinearly_constrained_optimization) and can be cited with the following BibTeX entry:

```
@unpublished{VanaretLeyffer2026,
  author = {Vanaret, Charlie and Leyffer, Sven},
  title = {Implementing a unified solver for nonlinearly constrained optimization},
  year = {2026},
  note = {Accepted to Mathematical Programming Computation on Feb 22, 2026}
}
```

## Credits

The theoretical abstract framework for unifying nonlinearly constrained optimization was developed by [Charlie Vanaret](https://github.com/cvanaret/) (Argonne National Laboratory & Zuse-Institut Berlin) and [Sven Leyffer](https://wiki.mcs.anl.gov/leyffer/index.php/Sven_Leyffer) (Argonne National Laboratory).
The interfaces and continuous integration infrastructure for Uno were developed and are maintained by [Alexis Montoison](https://github.com/amontoison) (Argonne National Laboratory) and Charlie Vanaret.
Uno itself was designed and implemented by Charlie Vanaret.
It is released under the MIT license (see the [license file](LICENSE)).

The contributors are (in alphabetical order):
[Oscar Dowson](https://github.com/odow), [Marcel Jacobse](https://github.com/mjacobse), [Arnav Kapoor](https://github.com/arnavk23), [David Kiessling](https://github.com/david0oo), [Rujia Liu](https://github.com/rujialiu), [Stefano Lovato](https://github.com/stefphd), [Manuel Schaich](https://github.com/worc4021), [Silvio Traversaro](https://github.com/traversaro).

The Uno logo was created by Charlie Vanaret based on a [saddle point icon by luimonts](https://thenounproject.com/icon/saddle-point-258207/) (CC BY 3.0).
