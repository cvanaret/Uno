SHELL := /bin/bash

# specify local compilers and compiler flags
CC     = gcc
CPP    = g++ -std=c++17
F77    = gfortran

include Makefile.local

FFLAGS = $(PROFILER)
LFLAGS = $(PROFILER)

# local libraries
L = -pthread -ldl -lf2c -fopenmp -lm
# LIBCOIN   = -lOsiClp -lClp -lCoinUtils -lOsi
LIBLINALG = -lblas -llapack

INCLUDES = -I base -I base/mechanism -I base/strategy -I base/subproblem -I logger -I solvers/QP -I solvers/linear -I solvers/NLP
# -I $(BFGS_DIR)

EXEC_NAME = argonot

VPATH = base solvers interfaces

SRC_FILES = base/Argonot.cpp base/Iterate.cpp \
	base/Problem.cpp base/Utils.cpp base/Statistics.cpp \
	base/Matrix.cpp base/HessianEvaluation.cpp base/Constraint.cpp \
	base/strategy/GlobalizationStrategy.cpp base/strategy/GlobalizationStrategyFactory.cpp \
	base/strategy/Filter.cpp base/strategy/FilterStrategy.cpp base/strategy/PenaltyMeritFunction.cpp \
	base/mechanism/GlobalizationMechanism.cpp base/mechanism/GlobalizationMechanismFactory.cpp \
	base/mechanism/TrustRegion.cpp base/mechanism/LineSearch.cpp \
	base/subproblem/Subproblem.cpp base/subproblem/SubproblemFactory.cpp base/subproblem/Direction.cpp \
	base/subproblem/ActiveSetMethod.cpp base/subproblem/SQP.cpp base/subproblem/SLP.cpp base/subproblem/Sl1QP.cpp base/subproblem/SLPEQP.cpp base/subproblem/InteriorPoint.cpp \
	solvers/QP/QPSolverFactory.cpp solvers/linear/LinearSolverFactory.cpp

# base/mechanism/TrustLineSearch.cpp

######################
# External libraries #
######################

DEFINES =
LIBS = -lf2c

# AMPL library
ifneq ($(HAS_AMPL),)
SRC_FILES := $(SRC_FILES) interfaces/AMPL/AMPLModel.cpp
DEFINES := $(DEFINES) $(HAS_AMPL)
LIBS := $(LIBS) $(AMPL_DIR)/amplsolver.a
INCLUDES := $(INCLUDES) -I interfaces/AMPL -I $(AMPL_DIR)
endif

# BQPD library
ifneq ($(HAS_BQPD),)
SRC_FILES := $(SRC_FILES) solvers/QP/BQPDSolver.cpp
DEFINES := $(DEFINES) $(HAS_BQPD)
LIBS := $(LIBS) $(BQPD_DIR)/libbqpd.a
endif

# MA57 library
ifneq ($(HAS_MA57),)
SRC_FILES := $(SRC_FILES) solvers/linear/MA57Solver.cpp
DEFINES := $(DEFINES) $(HAS_MA57)
LIBS := $(LIBS) -L$(MA57_DIR) -lma57 -l$(METIS)
endif

# PARDISO library
ifneq ($(HAS_PARDISO),)
SRC_FILES := $(SRC_FILES) solvers/linear/PardisoSolver.cpp
DEFINES := $(DEFINES) $(HAS_PARDISO)
LIBS := $(LIBS) -L$(PARDISO_DIR) -lpardiso600-GNU720-X86-64
endif

# BFGS library
ifneq ($(HAS_BFGS),)
SRC_FILES := $(SRC_FILES) solvers/NLP/LBFGSB.cpp
DEFINES := $(DEFINES) $(HAS_BFGS)
endif

OBJ_FILES = $(subst .cpp,.o,$(SRC_FILES)) solvers/QP/wdotd.o #$(BFGS_DIR)/lbfgsb.o $(BFGS_DIR)/timer.o $(BFGS_DIR)/blas.o $(BFGS_DIR)/linpack.o

# base/strategy/Tube.cpp base/strategy/TubeStrategy.cpp
# base/subproblem/AugmentedLagrangian.cpp \
#$(info HAS_BQPD is $(HAS_BQPD))

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
	rm -f $(EXEC_NAME) $(OBJ_FILES) main.o *~

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
