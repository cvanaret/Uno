# UNO (Unifying Framework for Optimization)

UNO is a C++ framework aiming at unifying most of the methods for solving nonlinearly constrained optimization problems.

## Table of contents
* [Installation instructions](#installation-instructions)
* [Compilation instructions](#compilation-instructions)
* [Unit tests](#unit-tests)
* [Autocompletion](#autocompletion)

## Installation instructions

* download optional interfaces:
  * AMPL: http://www.netlib.org/ampl/solvers/

* download optional solvers:
  * BQPD: https://www.mcs.anl.gov/~leyffer/solvers.html
  * MA57: http://www.hsl.rl.ac.uk/catalogue/ma57.html
  * PARDISO: https://www.pardiso-project.org/

* install BLAS, LAPACK and f2c:
```
sudo apt-get install libblas-dev liblapack-dev libf2c2-dev
```
* install cmake and ccmake:
```
sudo apt-get install cmake cmake-curses-gui
```

## Compilation instructions
1. Create a build directory in src/:
```
  cd src/ && mkdir build && cp uno.cfg build/
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
5. Compile in parallel (n being the number of threads, e.g. 6):
```
  make -jn
```
6. Run UNO:
```
  ./uno your_problem.nl
```

To compile the code with different configurations, simply create a build directory for each configuration and perform instructions 1 to 6.

## Unit tests
7. Install the GoogleTest suite:
```
  sudo apt-get install googletest
```
8. Perform steps 2 and 3
9. Run the test suite:
```
  ./run_unit_tests
```

## Autocompletion
To benefit from autocompletion, install the file uno-completion.bash:
```
  sudo cp uno-completion.bash /etc/bash_completion.d/
```
and open a new terminal.
