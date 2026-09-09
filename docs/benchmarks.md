# Optimization benchmarks

## Mittelmann benchmark

Uno v2.9.0 solves **45/47 instances** in the [Mittelmann benchmark](https://plato.asu.edu/ftp/ampl-nlp.html), in about the same time as IPOPT.

## CUTE benchmark

Uno presets have been tested against state-of-the-art solvers on 429 small problems of the CUTE benchmark [translated to AMPL](https://arnold-neumaier.at/glopt/coconut/Benchmark/Library2_new_v1.html).
The figure below (dated September 9, 2026) is a performance profile of Uno (solid lines) vs state-of-the-art solvers filterSQP and IPOPT (dotted lines); it shows how many problems are solved for a given budget of function evaluations (1 time, 2 times, 4 times, ..., $2^x$ times the number of objective evaluations of the best solver for each instance).

<p align="center">
   <img src="../figures/performance_profile_presets.png" alt="Performance profile of Uno vs filterSQP and IPOPT" width="75%" />
</p>

The performance profile below features the state-of-the-art solvers SNOPT, MINOS, LANCELOT, LOQO and CONOPT.

<p align="center">
   <img src="../figures/performance_profile_all.png" alt="Performance profile of Uno vs all solvers" width="75%" />
</p>

All log files can be found [here](https://github.com/cvanaret/nonlinear_optimization_solver_benchmark).