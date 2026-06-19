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
  - name: Charlie Vanaret^[corresponding author]
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
date: 8 June 2026
bibliography: paper.bib

---

# Summary

Uno is a composable software framework for nonlinearly constrained optimization written in modern C\texttt{++}.
It unifies the workflows of Lagrange-Newton methods, i.e., gradient-based algorithms that iteratively solve the KKT optimality conditions using Newton's method.
The central idea is to decompose these methods into reusable, interchangeable components (constraint reformulation, step computation, globalization techniques, and acceptance criteria) so that classical and hybrid algorithms can be assembled, compared, and tested within a single framework rather than reimplemented as separate solvers.
As of April 2026, Uno supports sequential (convex and nonconvex) quadratic programming, interior-point (barrier) methods, sequential linear programming, and unconstrained optimization.
For full mathematical details of the algorithms implemented in Uno, see [@VanaretLeyffer2026].

Uno has interfaces to Julia, Python, C, Fortran, and AMPL, enabling interoperability across scientific computing environments.
Precompiled artifacts are available on GitHub, and the solver can be accessed directly via `UnoSolver.jl` in Julia or `unopy` in Python.

# Statement of need

Nonlinearly constrained optimization is central to engineering, optimal control, machine learning, and scientific modeling [@nocedal2006], and plays a key role in mixed-integer nonlinear optimization [@lee2011].
The major solution paradigms (sequential quadratic programming, interior-point methods, and augmented Lagrangian methods) are usually developed and implemented as separate solver families.
Yet they share the same building blocks: step computation, constraint handling, Hessian models, and globalization strategies.
Because these building blocks are rigid inside monolithic codes, algorithmic ideas are typically tested at the level of *complete solvers* rather than *building blocks*: evaluating a new globalization strategy or Hessian approximation means reimplementing an entire solver or intrusively modifying an existing one.

Uno targets two audiences.
For **algorithm developers and optimization researchers**, it makes the building blocks explicit and recombinable across paradigms, so that a single new strategy can be swapped into an otherwise state-of-the-art method and benchmarked in isolation.
For **practitioners and end users**, it exposes those building blocks, pre-assembled into presets that reproduce established solvers, behind multiple language interfaces.
The gap Uno fills is the absence of a framework in which *production-grade* algorithmic components can be composed across solver paradigms under a unified abstraction.

# State of the field

Popular solvers such as IPOPT [@wachter2006implementation], KNITRO [@byrd2006], SNOPT [@gill2005], and WORHP [@buskens2012esa] are robust and efficient but  monolithic: parameter tuning is exposed while internal components remain rigid.
They are designed for end users rather than algorithmic experimentation, and testing a new variant usually requires intrusive modification or reimplementation.
A notable exception among production codes is MadNLP.jl [@shin2024accelerating], an actively developed Julia solver, though it remains focused on interior-point methods.

Frameworks that do emphasize modularity occupy a different niche.
General-purpose libraries such as `scipy.optimize`, NLopt, pyOpt, and pyOptSparse let users select among complete solvers, but they do not expose the *internal* algorithmic components for recombination.
Closest in spirit is modOpt [@joshy2026modopt], a Python environment that builds optimization algorithms from self-contained modules (line searches, Hessian approximations, merit functions) with an educational focus and transparent pedagogical implementations.

Uno is distinct on three counts.
First, its components are *production-grade* algorithmic building blocks compiled in modern C\texttt{++}, not pedagogical or wrapper modules.
Second, it composes them *across* paradigms -- SQP, interior-point, and SLP -- under one generic Lagrange-Newton abstraction, rather than within a single algorithm family.
Third, its presets reproduce state-of-the-art solvers at competitive performance (see below), so a component swapped into Uno is measured against a credible baseline rather than a teaching reference.

# Software design

Within Uno, strategies are organized into a coherent hierarchy, illustrated in the wheel of strategies (\autoref{fig:wheel}).
The outer ring lists the high-level *layers* of a Lagrange-Newton method (for example, constraint reformulation and globalization); the middle ring lists the algorithmic *ingredients* that Uno combines automatically to realize each layer; and the inner ring lists the concrete *strategies* available for each ingredient.
Strategies currently implemented in Uno are highlighted in green.

![Unification framework: wheel of strategies.\label{fig:wheel}](figures/wheel.pdf){ width=70% }

In Uno's object-oriented architecture, the ingredients of \autoref{fig:wheel} are abstract classes whose interfaces are implemented by subclasses (the strategies).
Uno's simplified UML diagram is shown in \autoref{fig:umldiagram}: inheritance is represented as dashed lines with white arrows, composition as solid lines with black diamonds; abstract classes are written in italic and subclasses in bold.
This architecture separates the mathematical logic of the algorithm (problem reformulation, subproblem definition, globalization) from the computational aspects (evaluating the model's functions, solving the subproblems).

![Uno's UML diagram.\label{fig:umldiagram}](figures/uml_diagram.pdf){ width=95% }

Uno implements a generic Lagrange-Newton method in which the abstract classes interact and exchange data while remaining agnostic to the underlying strategies.
Strategies are selected by the user at runtime via options.
Uno also provides *presets*: particular combinations of strategies and hyperparameters that correspond to state-of-the-art solvers.
Two presets are currently available: an \texttt{ipopt} preset mimicking the IPOPT solver [@wachter2006implementation] (a *line-search restoration filter interior-point method with exact Hessian and primal-dual inertia correction*), and a \texttt{filtersqp} preset mimicking the filterSQP solver [@fletcher1998user] (a *trust-region restoration filter SQP method with exact Hessian and no inertia correction*).
In principle all combinations of strategies can be generated, although some are not yet supported (e.g., an interior-point method with a trust-region constraint).

These two presets perform on a par with the state-of-the-art solvers IPOPT (Uno is slightly less robust) and filterSQP (Uno is slightly more robust) in terms of function evaluations on a set of 429 small CUTE instances with fewer than 100 variables and constraints [@bongartz1995cute] converted to AMPL.
An up-to-date performance profile is maintained on Uno's documentation page.
This benchmark shows that Uno's composability does not come at the cost of performance: assembling a method from modular components reproduces the behavior of the hand-written solvers it mimics.

Subproblem solvers are treated as interchangeable components that can be plugged in or swapped out without modifying the algorithmic logic, letting users match the solver to their problem structure, licensing constraints, or performance requirements.
Uno currently interfaces several established LP, QP, and linear solvers: BQPD [@fletcher2000stable], HiGHS [@huangfu2018parallelizing], MUMPS [@amestoy2000mumps], MA27 [@duffma27], MA57 [@duff2004ma57], SSIDS [@hogg2016sparse], and `Krylov.jl` [@montoison2023krylov].

# Interfaces

To make Uno accessible across scientific computing environments, we provide five language interfaces, each oriented toward a different user workflow.
The C interface is the foundation: it expresses the optimization process as the interaction of two opaque objects: a *model* (variables, bounds, constraints, objective, and derivative information) and a *solver* (algorithmic configuration, callbacks, and result such as primal-dual solutions, residuals, and iteration counts), and underpins the other interfaces.

| Interface                | Distribution                              | Typical workflow                                                       |
|--------------------------|-------------------------------------------|------------------------------------------------------------------------|
| AMPL [@gay1997hooking]   | standalone binary reading `.nl` files     | solve AMPL [@fourer1990] models without touching the C\texttt{++} core |
| C                        | core API (model + solver objects)         | embed Uno in C/C\texttt{++} applications; basis for other bindings     |
| Fortran                  | source bindings via `iso_c_binding`       | call Uno from Fortran codes across compilers and platforms             |
| Julia (`UnoSolver.jl`)   | registered package, precompiled artifacts | solve `JuMP.jl` and `NLPModels.jl` models, plug-and-play               |
| Python (`unopy`)         | PyPI wheels (`pip install unopy`)         | scripting and Python optimization workflows                            |

The Julia interface integrates Uno into the Julia ecosystem through a thin wrapper over the C API, an `NLPModels.jl` interface (e.g., `ADNLPModels.jl`, `ExaModels.jl`), and a `MathOptInterface.jl` interface for JuMP models.
A MATLAB interface is also under development.
Looking ahead, we plan to integrate Uno as a continuous relaxation solver within MINLP frameworks such as SCIP [@bolusani2024scip], which will require efficient reoptimization that reuses internal solver state (active sets, factorizations, preallocated workspace) across structurally similar solves.

# Research impact statement

Uno has demonstrated value both as a research tool and as a dependency of established software.
Its original authors used it to prototype a novel globalization strategy, the funnel method [@kiessling2025unified].
It has since been used as the lower-level solver in derivative-free bilevel optimization [@cesaroni2026derivative], been benchmarked on Euclidean problems with bound constraints [@baran2026riemannian], and been cited in [@gruss2024estimation; @desef2025optimization; @cederberg2026disciplined].

Uno is currently available as a nonlinear optimization solver in:

- Julia:
	* the `JuMP.jl` ecosystem,
	* `control-toolbox`, a collection of Julia packages for mathematical control and its applications,
- C++:
  * [CasADi](https://github.com/casadi/casadi), an open-source tool for nonlinear optimization and algorithmic differentiation [@Andersson2019],
- Python:
	* [pyOptSparse](https://github.com/mdolab/pyoptsparse), an object-oriented framework for formulating and solving nonlinear constrained optimization problems,
	* [DNLP](https://github.com/cvxgrp/DNLP), an extension of CVXPY to general nonlinear programming,
- Fortran:
	* [CUTEst](https://github.com/ralna/cuTEst), the Constrained and Unconstrained Testing Environment with safe threads for optimization software,
	* IMPL © /IMPL-DATA © by Industrial Algorithms Limited, a modeling and solving platform used in the process industries.

# Acknowledgments

This work was supported by the Applied Mathematics activity of the U.S. Department of Energy Office of Science, Advanced Scientific Computing Research, under Contract No. DE-AC02-06CH11357.
This work was also supported in part by NSF CSSI Grant No. 2104068.

# AI usage disclosure

No generative AI was used in the creation of the software, interfaces, implementation of algorithms, or documentation.
ChatGPT was used to check spelling, grammar, and clarity of the English text in this paper.

# References
