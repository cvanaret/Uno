# Uno, a unified solver for nonlinearly constrained optimization

**Uno** (Unifying Nonlinear Optimization)[^1] is a modular open-source solver for nonlinearly constrained optimization.
It unifies most derivative-based iterative (Lagrange-Newton) methods and organizes common algorithmic components (such as constraint reformulation, step computation, and globalization) into a coherent hierarchy:

<p align="center">
   <img src="figures/wheel.png" alt="Unifying framework: the wheel of strategies." width="60%" />
</p>

Uno allows the automatic generation of various strategy combinations on the fly with no programming effort from the user. While all combinations do not lead to convergent methods, some of them result in efficient solvers that may not exist as software implementations.

For a gentle introduction, watch the talk I gave at [JuMP-dev 2026](https://jump.dev/meetings/jumpdev2026/):
<div class="video-wrapper">
   <iframe src="https://www.youtube-nocookie.com/embed/hDPCVBQm8TE" title="Video title" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen>
   </iframe>
</div>

## Available methods

Uno currently implements the following strategies:

- **constraint relaxation strategies**: feasibility restoration;
- **inequality handling methods**: inequality constrained method, interior-point method;
- **Hessian models**: exact, L-BFGS, L-SR1, identity, zero;
- **inertia control strategies**: primal, primal-dual, none;
- **globalization strategies**: filter method, funnel method, merit function;
- **globalization mechanisms**: backtracking line search, trust-region method.

Uno also implements [presets](presets.md), that is strategy combinations that mimic existing solvers:

* `filtersqp` mimics filterSQP (trust-region feasibility restoration filter SQP method with exact Hessian);
* `ipopt` mimics IPOPT (line-search feasibility restoration filter barrier method with exact Hessian and primal-dual inertia correction).

The default preset `auto` decides between `filtersqp` and `ipopt`, depending on the properties of the problem.

## Who uses Uno?

Uno is currently available as a nonlinear optimization solver in:

- [JuMP.jl](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers)
- [CasADi](https://github.com/casadi/casadi)
- [pyOptSparse](https://github.com/mdolab/pyoptsparse)
- [DNLP](https://github.com/cvxgrp/DNLP), an extension of [CVXPY](https://www.cvxpy.org/) to general nonlinear programming
- [CUTEst](https://github.com/ralna/cuTEst/), the Constrained and Unconstrained Testing Environment with safe threads (CUTEst) for optimization software 
- [Vecchia.jl](https://github.com/cgeoga/Vecchia.jl), a package for Gaussian processes approximation
- [control-toolbox](https://github.com/control-toolbox), a collection of Julia packages for mathematical control and its applications
- [FelooPy](https://www.linkedin.com/posts/k-tafakkori_optimization101-operationsresearch-decisionscience-activity-7397646574035697664-AzmK), a user-friendly tool for coding, modeling, and solving decision problems
- [IMPL &copy; /IMPL-DATA &copy;](https://www.linkedin.com/posts/jeffrey-dean-kelly-a5420a6a_releases-cvanaretuno-activity-7388564004585160704-WSxz/) by [Industrial Algorithms Limited](https://www.industrialgorithms.ca/), a modeling and solving platform used in the process industries especially suited for economic, efficiency and emissions optimization and estimation

and more to come:

- [Pyomo](https://github.com/cvanaret/Uno/issues/319)
- [Minotaur](https://github.com/cvanaret/Uno/issues/107)
- [NEOS Server](https://neos-server.org/neos/solvers/)

[^1]: Uno was first introduced at the ISMP 2018 conference under the name Argonot.
