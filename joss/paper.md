---
title: 'UNO: A composable framework for nonlinear constrained optimization'
tags:
  - nonlinear programming
  - constrained optimization
  - SQP
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

UNO is a composable software framework for nonlinear constrained optimization written in modern C++.
It implements a unified Lagrange-Newton approach by decomposing classical algorithms into reusable components.
Supported methods include Sequential Quadratic Programming (SQP), interior-point barrier methods, and augmented Lagrangian approaches.

Optimization algorithms in UNO are structured as compositions of building blocks: search directions, constraint treatments, globalization strategies, and acceptance mechanisms.
This design allows classical and hybrid methods to be expressed, configured, and compared within the same framework.

The core C++ implementation separates mathematical concepts from implementation details.
Language bindings for Julia, Python, C, and Fortran make UNO accessible across scientific computing environments. Precompiled libraries and executables are provided on GitHub, and the solver can be used directly via `UnoSolver.jl` in Julia or `unopy` in Python.

# Statement of need

Nonlinear constrained optimization is central to engineering, optimal control, machine learning, and scientific modeling.
Existing solvers are typically monolithic, exposing only parameter tuning rather than structural composition.

This creates limitations for algorithmic research:

* Implementing new methods requires modifying complex legacy code.
* Comparisons of algorithmic strategies are difficult to reproduce.
* Hybrid methods combining multiple paradigms are hard to prototype.
* Teaching algorithmic components is challenging.

UNO addresses these gaps by enabling users to assemble algorithms from fundamental building blocks corresponding to mathematical concepts such as step computation, constraint modeling, penalty or barrier handling, and globalization strategies.
It serves both as a practical solver and as a platform for research and education.

# State of the field

Nonlinear constrained optimization solvers typically implement:

* Sequential Quadratic Programming
* Interior-point barrier methods
* Augmented Lagrangian methods

Existing solvers are robust but structurally rigid, which makes it difficult to exchange components or combine methods.
Modeling frameworks simplify problem formulation but rely on external solvers and do not expose internal algorithmic structure.
Research codes sometimes implement specialized methods, but components are not interchangeable, globalization strategies are hard-wired, and combining paradigms requires significant reimplementation.

UNO differs by expressing algorithms as compositions of reusable modules.
Classical solvers emerge as configurations within the same framework rather than independent implementations.
This enables reproducible research, controlled comparisons, and systematic exploration of hybrid methods.

# Software design

...

# Research impact statement

...

UNO is also a production-quality solver.
Its C++ core is designed to be efficient, and language bindings for Julia, Python, C, and Fortran make it directly usable in scientific and engineering workflows.
By representing classical algorithms within a single framework, experiments isolate algorithmic ideas from implementation details.
UNO functions as both a state-of-the-art solver and a platform for algorithmic research.

# References
