# Uno

<div align="center">

   *A modern, modular solver for nonlinearly constrained optimization*

</div>

Uno (Unifying Nonlinear Optimization) is a C++ library that unifies methods for solving nonlinearly constrained optimization problems of the form:

$$
\begin{align}
\min_{x \in \mathbb{R}^n}  & ~f(x) \\
\text{s.t.}                & ~c_L \le c(x) \le c_U \\
                           & ~x_L \le x \le x_U \\
\end{align}
$$

The theoretical abstract framework for unifying nonlinearly constrained optimization was developed by [Charlie Vanaret](https://github.com/cvanaret/) (Argonne National Laboratory & Zuse-Institut Berlin) and [Sven Leyffer](https://wiki.mcs.anl.gov/leyffer/index.php/Sven_Leyffer) (Argonne National Laboratory). Uno was designed and implemented by Charlie Vanaret. It is released under the MIT license (see the [license file](LICENSE)).

The contributors are (in alphabetical order):
Oscar Dowson [@odow](https://github.com/odow), David Kiessling [@david0oo](https://github.com/david0oo), Alexis Montoison [@amontoison](https://github.com/amontoison), Manuel Schaich [@worc4021](https://github.com/worc4021), Silvio Traversaro [@traversaro](https://github.com/traversaro), Rujia Liu [@rujialiu](https://github.com/rujialiu).

![Unit tests on Ubuntu workflow](https://github.com/cvanaret/Uno/actions/workflows/unit-tests-ubuntu.yml/badge.svg)

## Unifying nonlinearly constrained optimization

We argue that most optimization methods can be broken down into the following generic ingredients:
* a **constraint relaxation strategy**: a systematic way to relax the general constraints;
* an **inequality handling method**: a systematic way to handle the inequality constraints;
* a **Hessian model**: a model of the Lagrangian Hessian of the original problem;
* a **regularization strategy**: a strategy to regularize the Lagrangian Hessian or the augmented system of the reformulated problem;
* a **globalization strategy**: an acceptance test of the trial iterate;
* a **globalization mechanism**: a recourse action upon rejection of the trial iterate.

The following graph gives an overview of state-of-the-art strategies:
<p align="center">
   <img src="docs/figures/wheel.png" alt="Uno hypergraph" width="65%" />
</p>

**Any strategy combination** can be automatically generated without any programming effort from the user. Note that all combinations do not necessarily result in sensible algorithms, or even convergent approaches. For more details, check out our [preprint](https://www.researchgate.net/publication/381522383_Unifying_nonlinearly_constrained_nonconvex_optimization) or my [latest slides](https://www.researchgate.net/publication/390271091).

Uno implements **presets**, that is strategy combinations that correspond to existing solvers (as well as hyperparameter values found in their documentations):
* `filtersqp` mimics filterSQP (trust-region feasibility restoration filter SQP method with exact Hessian);
* `ipopt` mimics IPOPT (line-search feasibility restoration filter barrier method with exact Hessian and primal-dual regularization);
* `byrd` mimics Byrd's S $\ell_1$ QP (line-search $\ell_1$ merit S $\ell_1$ QP method with exact Hessian and primal regularization).

## Latest results (September 26, 2024)

Some of Uno combinations that correspond to existing solvers (called presets, see below) have been tested against state-of-the-art solvers on 429 small problems of the [CUTEst benchmark](https://arnold-neumaier.at/glopt/coconut/Benchmark/Library2_new_v1.html).
The figure below is a performance profile of Uno and state-of-the-art solvers filterSQP, IPOPT, SNOPT, MINOS, LANCELOT, LOQO and CONOPT; it shows how many problems are solved for a given budget of function evaluations (1 time, 2 times, 4 times, ..., $2^x$ times the number of objective evaluations of the best solver for each instance).

<p align="center">
   <img src="docs/figures/uno_performance_profile.png" alt="Performance profile of Uno" width="75%" />
</p>

All log files can be found [here](https://github.com/cvanaret/nonconvex_solver_comparison).

## How to cite Uno

### In an article

We have submitted our paper to the Mathematical Programming Computation journal. The preprint is available on [ResearchGate](https://www.researchgate.net/publication/381522383_Unifying_nonlinearly_constrained_nonconvex_optimization).

Until it is published, you can use the following bibtex entry:

```
@unpublished{VanaretLeyffer2024,
  author = {Vanaret, Charlie and Leyffer, Sven},
  title = {Unifying nonlinearly constrained nonconvex optimization},
  year = {2024},
  note = {Submitted to Mathematical Programming Computation}
}
```

### On social media

To mention Uno on Twitter, use [@UnoSolver](https://twitter.com/UnoSolver).  
To mention Uno on LinkedIn, use [#unosolver](https://www.linkedin.com/feed/hashtag/?keywords=unosolver).  

## Installation instructions

See the [INSTALL](INSTALL.md) file.

## Solving a problem with Uno

### Controlling Uno via options

Options can be set in three different ways (with decreasing precedence):
- passing an option file (`option_file=file`) that contains `option value` on each line;
- setting a preset that mimics an existing solver (`preset=[filtersqp|ipopt|byrd]`);
- setting individual options (see the [default options](https://github.com/cvanaret/Uno/blob/main/uno/options/DefaultOptions.cpp)).

### Interfaces

#### AMPL/nl files
To solve an AMPL model in the [.nl format](https://en.wikipedia.org/wiki/Nl_(format)), type in the `build` directory: ```./uno_ampl model.nl -AMPL [option=value ...]```  
where ```[option=value ...]``` is a list of options separated by spaces. 

A couple of CUTEst instances are available in the `/examples` directory.

#### Julia
Uno can be installed in Julia via [Uno_jll.jl](https://github.com/JuliaBinaryWrappers/Uno_jll.jl) and used via [AmplNLWriter.jl](https://juliahub.com/ui/Packages/General/AmplNLWriter.jl). An example can be found [here](https://discourse.julialang.org/t/the-uno-unifying-nonconvex-optimization-solver/115883/15?u=cvanaret).

### Combining strategies on the fly

For an overview of the available strategies, type: ```./uno_ampl --strategies```:
- to pick a constraint relaxation strategy, use the argument: ```constraint_relaxation_strategy=[feasibility_restoration|l1_relaxation]```  
- to pick an inequality handling method, use the argument: ```subproblem=[QP|LP|primal_dual_interior_point]``` 
- to pick a Hessian model, use the argument: ```hessian_model=[exact|identity|zero]``` 
- to pick a regularization strategy, use the argument: ```regularization_strategy=[primal|primal_dual|none]``` 
- to pick a globalization strategy, use the argument: ```globalization_strategy=[l1_merit|fletcher_filter_method|waechter_filter_method|funnel_method]```   
- to pick a globalization mechanism, use the argument : ```globalization_mechanism=[LS|TR]```  
