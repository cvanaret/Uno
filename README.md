# What is Uno?
Uno (Unifying Nonlinear Optimization) is a C++ framework aiming at unifying most of the methods for solving nonlinearly constrained optimization problems.

Check out my [presentation at the ICCOPT 2022 conference](https://www.researchgate.net/publication/362254109).
This is joint work with Sven Leyffer (Argonne National Laboratory).

# How to cite Uno

Please be patient, we are actively working on our article.

# Contributions

Developed by Charlie Vanaret (Technische Universität Berlin)

# License

Uno is released under the MIT license (see the [license file](LICENSE)).

# Installation instructions

* download the AMPL solver library (ASL): http://www.netlib.org/ampl/solvers/

* download optional solvers:
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

# Compilation instructions
1. Create a build directory in src/:
```
  cd src/ && mkdir build
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

To compile the code with different configurations, simply create a build directory for each configuration and perform instructions 1 to 5.

# Solving a problem with Uno
To solve an AMPL model, type:
```
./uno_ampl path_to_file/file.nl
```
To choose a globalization mechanism, use the argument:
```
-mechanism [LS|TR]
```
To choose a constraint relaxation strategy, use the argument:
```
-constraint-relaxation [feasibility-restoration|l1-relaxation]
```
To choose a globalization strategy, use the argument:
```
-strategy [penalty|filter|nonmonotone-filter]
```
To choose a subproblem method, use the argument:
```
-subproblem [QP|LP|barrier]
```
To choose a preset, use the argument:
```
-preset [filtersqp|ipopt|byrd]
```
The options can be combined in the same command line. Autocompletion is active.

# Unit tests
7. Install the GoogleTest suite:
```
  sudo apt-get install googletest
```
8. Perform steps 2 and 3
9. Run the test suite:
```
  ./run_unotest
```

## Autocompletion
To benefit from autocompletion, install the file uno_ampl-completion.bash:
```
  sudo cp uno_ampl-completion.bash /etc/bash_completion.d/
```
and open a new terminal.
