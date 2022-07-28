# Uno (Unifying Nonlinear Optimization)

## What is Uno?

Uno (Unifying Nonlinear Optimization) is a C++ library aiming at unifying most of the methods for solving nonlinearly constrained optimization problems of the form:

```
   min     f(x)
  x ∈ Rⁿ

   s.t.    c_L ≤ c(x) ≤ c_U
           x_L ≤  x   ≤ x_U
```

Uno implements an abstract framework based on four ingredients:
* **constraint relaxation strategy**: a systematic way to relax the nonlinear constraints;
* **subproblem**: a local model of the (possibly relaxed) problem at the current primal-dual iterate;
* **globalization strategy**: an acceptance test of the trial iterate;
* **globalization mechanism**: a recourse action upon rejection of the trial iterate.

The following hypergraph illustrates how state-of-the-art solvers can be decomposed in terms of the four ingredients.
<p align="center">
   <img src="docs/figures/combination_hypergraph.png" alt="Combination hypergraph" width="75%" />
</p>

Uno 1.0 implements the following strategies. Any strategy combination can be generated without any programming effort from the user. Note that all combinations do not necessarily result in sensible algorithms, or even convergent approaches.
<p align="center">
   <img src="docs/figures/hypergraph_uno.png" alt="Uno 1.0 hypergraph" width="65%" />
</p>

Check out my [presentation at the ICCOPT 2022 conference](https://www.researchgate.net/publication/362254109).
This is joint work with Sven Leyffer (Argonne National Laboratory).

## How to cite Uno

Please be patient, we are actively working on our article.

## Contributions

Developed by Charlie Vanaret (Technische Universität Berlin)

## License

Uno is released under the MIT license (see the [license file](LICENSE)).

## Installation instructions

### Packages and libraries

* download the AMPL solver library (ASL): http://www.netlib.org/ampl/solvers/

* download **optional** solvers:
    * BQPD: https://www.mcs.anl.gov/~leyffer/solvers.html
    * MA57: http://www.hsl.rl.ac.uk/catalogue/ma57.html

* install BLAS, LAPACK and f2c:
```
sudo apt-get install libblas-dev liblapack-dev libf2c2-dev
```
* install cmake and ccmake (CMake curses interface):
```
sudo apt-get install cmake cmake-curses-gui
```

### Compilation

1. Create a `build` directory in the main directory:
```
  mkdir build
```
2. Move to the build directory:
```
  cd build/
```
3. Type cmake:
```
  cmake ..
```
4. Use ccmake to provide the paths to the required and optional libraries:
```
  ccmake ..
```
5. Compile in parallel (`n` being the number of threads, e.g. 6):
```
  make -jn
```
6. To print the version, type:
```
  ./uno_ampl -v
```

To compile the code with different configurations, simply create a `build` directory for each configuration and perform instructions 1 to 5.

### Unit tests

7. Install the GoogleTest suite:
```
  sudo apt-get install googletest
```
8. Perform steps 2 and 3
9. Run the test suite:
```
  ./run_unotest
```

### Autocompletion

To benefit from autocompletion, install the file `uno_ampl-completion.bash`:
```
  sudo cp uno_ampl-completion.bash /etc/bash_completion.d/
```
and open a new terminal.

## Solving a problem with Uno

To solve an AMPL model, type in the `build` directory:
```
./uno_ampl path_to_file/file.nl
```

### Combination of ingredients

To pick a globalization mechanism, use the argument (choose one of the possible options in brackets):
```
-mechanism [LS|TR]
```
To pick a constraint relaxation strategy, use the argument:
```
-constraint-relaxation [feasibility-restoration|l1-relaxation]
```
To pick a globalization strategy, use the argument:
```
-strategy [penalty|filter|nonmonotone-filter]
```
To pick a subproblem method, use the argument:
```
-subproblem [QP|LP|barrier]
```
The options can be combined in the same command line.

### Presets

Uno presets are strategy combinations that correspond to existing solvers (as well as known values for their hyperparameters). Uno 1.0 implements three presets:
* `filtersqp` mimics filterSQP (trust-region feasibility restoration filter SQP method);
* `ipopt` mimics IPOPT (line-search feasibility restoration filter barrier method);
* `byrd` mimics Byrd's S$\ell_1$QP (line-search $\ell_1$ merit S$\ell_1$QP method).

To pick a preset, use the argument:
```
-preset [filtersqp|ipopt|byrd]
```
