# Uno

<div align="center">

   *A modern, modular solver for nonlinearly constrained nonconvex optimization*

</div>

Uno (Unifying Nonlinear Optimization) is a C++ library that unifies methods for solving nonlinearly constrained optimization problems of the form:

$$
\begin{align}
\min_{x \in \mathbb{R}^n}  & ~f(x) \\
\text{s.t.}                & ~c_L \le c(x) \le c_U \\
                           & ~x_L \le x \le x_U \\
\end{align}
$$

The theoretical abstract framework for unifying nonlinearly constrained nonconvex optimization was developed by [Charlie Vanaret](https://github.com/cvanaret/) (Zuse-Institut Berlin) and [Sven Leyffer](https://wiki.mcs.anl.gov/leyffer/index.php/Sven_Leyffer) (Argonne National Laboratory).  

Uno was designed and implemented by Charlie Vanaret. [Silvio Traversaro](https://github.com/traversaro) contributed to the CMakeLists.

Uno is released under the MIT license (see the [license file](LICENSE)).

## Unifying nonlinearly constrained nonconvex optimization

We argue that most optimization methods can be broken down into four generic ingredients:
* a **constraint relaxation strategy**: a systematic way to relax the nonlinear constraints;
* a **subproblem**: a local model of the (possibly relaxed) problem at the current primal-dual iterate;
* a **globalization strategy**: an acceptance test of the trial iterate;
* a **globalization mechanism**: a recourse action upon rejection of the trial iterate.

<!--
The following hypergraph illustrates how some of the state-of-the-art solvers can be decomposed in terms of the four ingredients:
<p align="center">
   <img src="docs/figures/combination_hypergraph.png" alt="Combination hypergraph" width="75%" />
</p>

## Uno 1.0
-->

Uno 1.0 implements the following strategies:
<p align="center">
   <img src="docs/figures/hypergraph_uno.png" alt="Uno 1.0 hypergraph" width="65%" />
</p>

**Any strategy combination** can be automatically generated without any programming effort from the user. Note that all combinations do not necessarily result in sensible algorithms, or even convergent approaches. For more details, check out our [preprint](https://www.researchgate.net/publication/381522383_Unifying_nonlinearly_constrained_nonconvex_optimization) or my [presentation at the ICCOPT 2022 conference](https://www.researchgate.net/publication/362254109).

Uno 1.0 implements three **presets**, that is strategy combinations that correspond to existing solvers (as well as hyperparameter values found in their documentations):
* `filtersqp` mimics filterSQP (trust-region feasibility restoration filter SQP method);
* `ipopt` mimics IPOPT (line-search feasibility restoration filter barrier method);
* `byrd` mimics Byrd's S $\ell_1$ QP (line-search $\ell_1$ merit S $\ell_1$ QP method).

## Latest results (April 27, 2023)

Some of Uno combinations that correspond to existing solvers (called presets, see below) have been tested against state-of-the-art solvers on 429 small problems of the [CUTEst benchmark](https://arnold-neumaier.at/glopt/coconut/Benchmark/Library2_new_v1.html).
The figure below is a performance profile of Uno and state-of-the-art solvers filterSQP, IPOPT, SNOPT, MINOS, LANCELOT, LOQO and CONOPT; it shows how many problems are solved for a given budget of function evaluations (1 time, 2 times, 4 times, ..., $2^x$ times the number of objective evaluations of the best solver for each instance).

<p align="center">
   <img src="docs/figures/uno_performance_profile.png" alt="Performance profile of Uno 1.0" width="75%" />
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

At the moment, Uno only reads models from [.nl files](https://en.wikipedia.org/wiki/Nl_(format)). A couple of CUTEst instances are available in the `/examples` directory.

To solve an AMPL model, type in the `build` directory: ```./uno_ampl path_to_file/file.nl```

To use Uno with Julia/JuMP, a solution in the short term is to use the package [AmplNLWriter.jl](https://juliahub.com/ui/Packages/General/AmplNLWriter.jl) to dump JuMP models into .nl files.

### Combination of ingredients

To pick a globalization mechanism, use the argument (choose one of the possible options in brackets): ```-globalization_mechanism [LS|TR]```  
To pick a constraint relaxation strategy, use the argument: ```-constraint_relaxation_strategy [feasibility_restoration|l1_relaxation]```  
To pick a globalization strategy, use the argument: ```-globalization_strategy [l1_merit|fletcher_filter_strategy|waechter_filter_strategy]```  
To pick a subproblem method, use the argument: ```-subproblem [QP|LP|primal_dual_interior_point]```  
The options can be combined in the same command line.

For an overview of the available strategies, type: ```./uno_ampl --strategies```

To pick a preset, use the argument: ```-preset [filtersqp|ipopt|byrd]```
