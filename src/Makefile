SHELL := /bin/bash

# specify local compilers and compiler flags
PROFILER = #-pg
CC     = gcc
CPP    = g++ -std=c++11
F77    = gfortran
FFLAGS = $(PROFILER)

include Makefile.local

LFLAGS = $(PROFILER)

# local libraries
L = -pthread -ldl -lf2c
LIBCOIN   = -lOsiClp -lClp -lCoinUtils -lOsi
LIBLINALG = -lblas -llapack

INCLUDES = -I base -I base/mechanism -I base/strategy -I base/subproblem -I logger -I interfaces/AMPL -I solvers/QP -I solvers/KKT -I solvers/NLP -I $(AMPL_DIR)
# -I $(BFGS_DIR)

SRC_FILES = base/Argonot.cpp base/Iterate.cpp \
	base/Problem.cpp interfaces/AMPL/AMPLModel.cpp \
	base/Matrix.cpp base/HessianEvaluation.cpp base/Constraint.cpp \
	base/strategy/GlobalizationStrategy.cpp base/strategy/Filter.cpp base/strategy/FilterStrategy.cpp base/strategy/PenaltyMeritFunction.cpp \
	base/strategy/GlobalizationStrategyFactory.cpp \
	base/mechanism/GlobalizationMechanism.cpp base/mechanism/TrustRegion.cpp \
	base/mechanism/LineSearch.cpp base/mechanism/TrustLineSearch.cpp \
	base/mechanism/GlobalizationMechanismFactory.cpp \
	base/subproblem/Subproblem.cpp base/subproblem/ActiveSetMethod.cpp base/subproblem/SQP.cpp base/subproblem/SLP.cpp base/subproblem/Sl1QP.cpp base/subproblem/InteriorPoint.cpp \
	base/subproblem/SubproblemFactory.cpp base/subproblem/SubproblemSolution.cpp base/Utils.cpp \
	solvers/QP/BQPDSolver.cpp solvers/QP/QPSolverFactory.cpp \
	solvers/KKT/MA57Solver.cpp solvers/NLP/LBFGSB.cpp
OBJ_FILES = $(subst .cpp,.o,$(SRC_FILES)) solvers/QP/wdotd.o $(BFGS_DIR)/lbfgsb.o $(BFGS_DIR)/timer.o $(BFGS_DIR)/blas.o $(BFGS_DIR)/linpack.o
EXEC_NAME = argonot

#   base/subproblem/SLPEQP.cpp
# base/strategy/Tube.cpp base/strategy/TubeStrategy.cpp
# base/subproblem/AugmentedLagrangian.cpp \

DEFINES =# -D BQPD # https://stackoverflow.com/questions/10337248/what-does-the-g-d-flag-do
LIBS = $(AMPL_DIR)/amplsolver.a $(BQPD_DIR)/libbqpd.a -L$(MA57_DIR) -lma57 -l$(METIS) -lf2c

#%.o: %.cpp %.hpp
.cpp.o:
	$(CPP) $(CFLAGS) $(DEFINES) $(INCLUDES) -c $< -o $@ ;

.c.o:
	$(CC) -c $(FFLAGS) $(CFLAGS) $(DEFINES) $(INCLUDES) $*.c -o $@

.f.o:
	$(F77) -c $(FFLAGS) $< -o $@

.f90.o:
	$(F77) -I $(FFLAGS) -c $< -o $@ ;

all: $(EXEC_NAME) $(AMPL_DIR)/amplsolver.a

$(EXEC_NAME): main.o $(OBJ_FILES) $(AMPL_DIR)/amplsolver.a
	$(CPP) $(LFLAGS) -o $(EXEC_NAME) main.o $(OBJ_FILES) $(LIBS) $(LIBLINALG) $(L) $(LIBCOMPIL)

$(AMPL_DIR)/amplsolver.a:
	cd $(AMPL_DIR); make amplsolver.a

clean: cleanhere

cleanhere:
	rm -f $(EXEC_NAME) $(OBJ_FILES) *~

xsum.out: $SRC_FILES
	xsum $SRC_FILES >xsum1.out
	cmp xsum0.out xsum1.out && mv xsum1.out xsum.out || diff xsum[01].out

filter-al: filter-al.o $(OBJ_FILES) $(AMPL_DIR)/amplsolver.a
	$(CPP) $(LFLAGS) -o filter-al filter-al.o $(OBJ_FILES) $(LIBS) $(LIBLINALG) $(L) $(LIBCOMPIL)

#depend: .depend

#.depend: $(SOURCES)
#	rm -f ./.depend
#	$(CPP) $(CFLAGS) -MM $^>>./.depend;

#include .depend
