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
  - name: Charlie Vanaret
    orcid: 0000-0002-1131-7631
    affiliation: 1
  - name: Alexis Montoison
    orcid: 0000-0002-3403-5450
    affiliation: 2
affiliations:
 - name: Applied Optimization Department, Zuse-Institut Berlin, Germany
   index: 1
 - name: Mathematics and Computer Science Division, Argonne National Laboratory, USA
   index: 2
date: 20 February 2026
bibliography: paper.bib
---

# Summary

Uno is a composable software framework for nonlinearly constrained optimization written in modern C\texttt{++}.
It unifies the workflows of Lagrange-Newton methods, i.e., gradient-based algorithms that iteratively solve the KKT optimality conditions using Newton's method.
As of February 2026, Uno supports sequential (convex and nonconvex) quadratic programming, interior-point (barrier) methods, and sequential linear programming.

Uno breaks down optimization algorithms into reusable modular components such as step computation, constraint reformulation, globalization techniques, and acceptance criteria.
This allows classical and hybrid methods to be configured and compared within a single framework.

The core C\texttt{++} code of Uno is organized into modular, object-oriented components that separate the mathematical logic of the algorithms from implementation details such as memory management, data structures, and computational routines.
For full mathematical details of the algorithms implemented in Uno, see [@VanaretLeyffer2024].
Uno provides interfaces to Julia, Python, C, Fortran, and AMPL, enabling use across scientific computing environments.
Precompiled artifacts are available on GitHub, and the solver is directly available via `UnoSolver.jl` in Julia or `unopy` in Python.

# Statement of need

Nonlinearly constrained optimization is central to engineering, optimal control, machine learning, and scientific modeling [@nocedal2006].
It also plays a key role in mixed-integer nonlinear optimization [@lee2011].

Popular solvers such as IPOPT [@wachter2006implementation], KNITRO [@byrd2006], and SNOPT [@gill2005] have proven robust and efficient, but they are typically monolithic: they expose parameter tuning while keeping internal algorithmic components rigid and inaccessible.
This limits algorithmic research, making it hard to prototype hybrid methods, to evaluate different approaches, or to teach the underlying techniques.

Uno addresses these gaps by providing a composable framework in which algorithms emerge from modular code components corresponding to mathematical concepts such as step computation, constraint reformulation, and globalization techniques.
It enables rapid prototyping of new methods and serves both research and educational purposes.

Typical nonlinear solvers implement strategies such as sequential quadratic programming, interior-point methods, and augmented Lagrangian methods.
In Uno's unification framework, these strategies are organized into a coherent hierarchy, as illustrated in the wheel of strategies (\autoref{fig:wheel}): the outer ring represents high-level layers, the middle ring represents algorithmic ingredients that are automatically combined within Uno, and the inner ring lists possible strategies for each of the ingredients.
The strategies currently implemented in Uno are highlighted in green.

![Unification framework: wheel of strategies.\label{fig:wheel}](figures/wheel.pdf){ width=70% }

# Software design

The architecture of Uno follows an object-oriented design in which the ingredients of \autoref{fig:wheel} are abstract classes whose interfaces should be implemented by subclasses (the strategies).
For instance, `BacktrackingLineSearch` and `TrustRegionMethod` both inherit from the abstract class `GlobalizationMechanism` and implement its interface.
Uno's simplified UML diagram is shown in \autoref{fig:umldiagram}.
Inheritance is represented as dashed lines with white arrows, while composition is represented as solid lines with black diamonds.
Abstract classes are written in italic, while subclasses are written in bold.

Uno implements a generic Lagrange-Newton method in which the abstract classes interact with one another and exchange data, while being agnostic of the underlying strategies.
Strategies are picked by the user at runtime via options.
This modular architecture offers a clear separation between the mathematical logic of the optimization algorithm (the reformulation of the problem, the definition of the subproblem, and the globalization techniques) and the computational aspects (evaluating the model's functions, and solving the subproblems).

Uno also implements presets, that is particular combinations of strategies (and sets of hyperparameters) that correspond to state-of-the-art solvers.
Uno currently implements an \texttt{ipopt} preset that mimics the IPOPT solver [@wachter2006implementation] (a *line-search restoration filter interior-point method with exact Hessian and primal-dual inertia correction*), and a \texttt{filtersqp} preset that mimics the filterSQP solver [@fletcher1998user] (a *trust-region restoration filter SQP method with exact Hessian and no inertia correction*).
While, in theory, all combinations of strategies may be generated, some are not supported yet (e.g., interior-point method with a trust-region constraint).

![Uno's UML diagram.\label{fig:umldiagram}](figures/uml_diagram.pdf){ width=95% } 

# Interfaces

To make Uno accessible to a wide range of users, we provide multiple language interfaces.

The AMPL Solver Library [@gay1997hooking] interface gives access to the AMPL [@fourer1990] modeling language for optimization.
It is distributed as a binary that takes a compiled AMPL model (.nl file) as input, allowing users to solve problems without directly interacting with the C\texttt{++} core.

The C interface provides direct access to Uno's core functionality while maximizing interoperability with other programming languages and tools.
The optimization process is expressed as the interaction between two independent opaque objects: a model describing the optimization problem and a solver holding the algorithmic configuration and execution state.
The interface is therefore centered around two main structures:

* **Model**: represents an optimization problem and stores information about variables, bounds, constraints, the objective function, and derivative information. Users create a model and set its components (objective, constraints, derivatives, initial primal-dual point).
* **Solver**: represents the algorithm used to solve a given model. Users set solver options, attach callbacks, and access results such as primal and dual solutions, residuals, iteration counts, and performance metrics.

The Fortran interface provides access to Uno through the C API using `iso_c_binding`.
It is split into two files: `uno_c.f90` for low-level C bindings, and `uno_fortran.f90` for Fortran-friendly wrappers that handle string arguments.
The interface can be included directly in a source file or wrapped in a module for cleaner `use` statements.
It is provided as source rather than a precompiled Fortran module (.mod) to maximize portability and interoperability across compilers and platforms.

The Julia interface is distributed as the registered package `UnoSolver.jl`.
It integrates Uno into the Julia optimization ecosystem through:

* a thin wrapper around the full C API,
* an interface to `NLPModels.jl` for solving problems following the NLPModels API, such as `ADNLPModels.jl` or `ExaModels.jl`,
* an interface to `MathOptInterface.jl` for handling JuMP models.

Precompiled artifacts are downloaded automatically, making the package plug-and-play without requiring user compilation.

The Python interface is registered on PyPI as `unopy`, providing Uno bindings via precompiled wheels on most platforms.

A MATLAB interface is also under development, further expanding Uno's accessibility.

# Research impact

Uno's \texttt{ipopt} and \texttt{filtersqp} presets currently perform on a par with the state-of-the-art solvers IPOPT (Uno is slightly less robust) and filterSQP (Uno is slightly more robust) in terms of function evaluations on a set of 429 small CUTE instances [@bongartz1995cute].
An up-to-date performance profile is maintained on Uno's GitHub `README` page.

Uno's ongoing developments were presented at several international conferences over the years (ISMP 2018, SIAM OP 2021, ICCOPT 2022, ISMP 2024, and ICCOPT 2025), and were the subject of invited talks at Zuse-Institut Berlin (2020), Argonne National Laboratory (2022), and KU Leuven (2025).
This resulted in scientific cooperations with the MECO team at KU Leuven and the HiGHS team in Edinburgh.

Uno is currently used as a nonlinear optimization solver in:

* `JuMP.jl` ecosystem,
* DNLP, an extension of CVXPY to general nonlinear programming,
* `Vecchia.jl`, a package for Gaussian processes approximation,
* FelooPy, a user-friendly tool for coding, modeling, and solving decision problems,
* IMPL © /IMPL-DATA © by Industrial Algorithms Limited, a modeling and solving platform used in the process industries especially suited for economic, efficiency and emissions optimization and estimation.

Ongoing discussions and community interest indicate potential future integrations in CasADi, Pyomo, pyOptSparse, Minotaur, and the NEOS Server.

# Acknowledgments

This work was supported by the Applied Mathematics activity of the U.S. Department of Energy Office of Science, Advanced Scientific Computing Research, under Contract No. DE-AC02-06CH11357.
This work was also supported in part by NSF CSSI Grant No. 2104068.

# AI usage disclosure

No generative AI was used in the creation of the software, interfaces, implementation of algorithms, or documentation.
ChatGPT was used to check spelling, grammar, and clarity of the English text in this paper.

# References
