## Installation instructions

### Packages and libraries

* download the AMPL solver library (ASL): http://www.netlib.org/ampl/solvers/

* download **optional** solvers:
    * BQPD (indefinite null-space QP solver): https://www.mcs.anl.gov/~leyffer/solvers.html
    * MA57 (sparse indefinite symmetric linear solver): http://www.hsl.rl.ac.uk/catalogue/ma57.html

* install BLAS and LAPACK: ```sudo apt-get install libblas-dev liblapack-dev```
* install cmake (and optionally ccmake, CMake curses interface): ```sudo apt-get install cmake cmake-curses-gui```

### Compilation

1. Create a `build` directory in the main directory: ```mkdir build```
2. Move to the build directory: ```cd build/```
3. Execute cmake by providing the paths to the required and optional libraries ASL, BQPD and MA57:  
```cmake -Dbqpd=path -Dma57=path -Damplsolver=path -DCMAKE_BUILD_TYPE=[Release|Debug] ..```
4. **(or)** Use ccmake to provide the paths to the required and optional libraries: ```ccmake ..```
5. Compile (in parallel: `n` being the number of threads, e.g. 6): ```make -jn```

To compile the code with different configurations, simply create a `build` directory for each configuration and perform instructions 1 to 5.

### Unit tests

6. Install the GoogleTest suite: ```sudo apt-get install googletest```
7. Perform steps 2 and 3 with the flag ```-DWITH_GTEST=ON```
8. Run the test suite: ```./run_unotest```

### Autocompletion

To benefit from autocompletion, install the file `uno_ampl-completion.bash`: ```sudo cp uno_ampl-completion.bash /etc/bash_completion.d/```
and open a new terminal.
