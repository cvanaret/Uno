---
title: 'Uno: a composable framework for nonlinearly constrained optimization'
tags:
  - nonlinear programming
  - constrained optimization
  - sequential quadratic programming
  - interior-point methods
  - Lagrange-Newton methods
  - optimization frameworks
  - scientific computing
  - Julia
  - Python
  - C
  - Fortran
  - AMPL
authors:
  - name: Alexis Montoison
    orcid: 0000-0002-3403-5450
    affiliation: 1
  - name: Charlie Vanaret
    orcid: 0000-0002-1131-7631
    affiliation: 2
affiliations:
 - name: Mathematics and Computer Science Division, Argonne National Laboratory, USA
   index: 1
 - name: Applied Optimization Department, Zuse-Institut Berlin, Germany
   index: 2
date: 13 February 2026
bibliography: paper.bib
---

# Summary

Uno is a composable software framework for nonlinearly constrained optimization written in modern C\texttt{++}. It unifies the workflows of Lagrange-Newton methods, i.e., gradient-based algorithms that iteratively solve the KKT optimality conditions with Newton's method. As of February 2026, Uno supports sequential (convex and nonconvex) quadratic programming, interior-point (barrier) methods, and sequential linear programming.

Uno breaks down optimization algorithms into reusable modular components such as step computation, constraint reformulation, globalization techniques, and acceptance criteria.
This allows classical and hybrid methods to be configured and compared within a single framework.

The core C\texttt{++} implementation separates mathematical abstractions from implementation details through a clear software abstraction layer.
Bindings for Julia, Python, C, and Fortran support the use of Uno across scientific computing environments.
Precompiled artifacts are available on GitHub.
The solver can be accessed directly via `UnoSolver.jl` in Julia or `unopy` in Python.

# Statement of need

Nonlinearly constrained optimization is central to engineering, optimal control, machine learning, and scientific modeling. It also plays a central role in mixed-integer nonlinear optimization (MINLP) in which a sequence of continuous relaxations is solved.
Existing nonlinear programming solvers are typically monolithic, exposing only parameter tuning rather than structural composition.

This creates limitations for algorithmic research:

* Implementing new methods requires modifying complex legacy code.
* The comparison of algorithmic strategies is difficult to reproduce.
* Hybrid methods combining multiple paradigms are hard to prototype.
* Teaching algorithmic components is challenging.

Uno addresses these gaps by enabling users to assemble algorithms from fundamental building blocks that correspond to mathematical concepts such as step computation, constraint reformulation (penalty or barrier), and globalization techniques.
Uno serves both as a practical solver and as a platform for research and education.

# State of the field

Nonlinearly constrained optimization solvers typically implement:

* Sequential Quadratic Programming
* Interior-point/barrier methods
* Augmented Lagrangian methods

Existing solvers are robust but structurally rigid, which makes it difficult to exchange components or combine methods.
Modeling frameworks simplify the problem formulation but rely on external solvers, and do not expose internal algorithmic structure.
Research codes sometimes implement specialized methods, but components are not interchangeable, globalization strategies are hard-wired, and combining paradigms requires significant reimplementation.

Uno differs in that algorithms emerge as possible combinations of strategies within a single composable framework, rather than independent implementations. This enables reproducible research, controlled comparisons, and systematic exploration of hybrid methods.
An extended wheel of strategies, organized as layers and ingredients, is shown in \autoref{fig:wheel}. Note that all strategies are not available in Uno yet. (TODO: use 2 different colors: one for implemented, one for not yet).

![Unification framework: wheel of strategies.\label{fig:wheel}](https://raw.githubusercontent.com/cvanaret/Uno/refs/heads/main/docs/figures/wheel.png){ width=60% }

# Software design

The architecture of Uno follows a usual object-oriented design in which abstract classes define interfaces that should be implemented by subclasses. For instance, `BacktrackingLineSearch` and `TrustRegionMethod` both inherit from the abstract class `GlobalizationMechanism` and implement its interface. 

# Research impact statement

...

Uno is also a production-quality solver.
Its C\texttt{++} core is designed to be efficient, and language bindings for Julia, Python, C, and Fortran make it directly usable in scientific and engineering workflows.
By representing classical algorithms within a single framework, experiments isolate algorithmic ideas from implementation details.
Uno functions as both a state-of-the-art solver and a platform for algorithmic research.

# References
