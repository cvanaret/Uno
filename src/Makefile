# specify local compilers and compiler flags
PROFILER =# -pg
CC     = gcc
CPP    = g++
F77    = gfortran
FFLAGS = $(PROFILER)

include Makefile.local

LFLAGS = $(PROFILER)

# local libraries
L = -ldl
LIBCOIN   = -lOsiClp -lClp -lCoinUtils -lOsi
LIBLINALG = #-lblas -llapack

INCLUDES = -I base  -I base/globalization_mechanism -I base/globalization_strategy -I base/local_approximation -I interfaces -I logger -I interfaces/AMPL -I solvers -I $(AMPL_DIR)

SRC_FILES = main.cpp base/Argonot.cpp base/Iterate.cpp \
	base/Problem.cpp interfaces/AMPL/AMPLModel.cpp \
	base/Matrix.cpp \
	base/globalization_strategy/GlobalizationStrategy.cpp base/globalization_strategy/TwoPhaseStrategy.cpp \
	base/globalization_strategy/Filter.cpp base/globalization_strategy/FilterStrategy.cpp \
	base/globalization_strategy/PenaltyStrategy.cpp \
	base/globalization_mechanism/GlobalizationMechanism.cpp base/globalization_mechanism/TrustRegion.cpp base/globalization_mechanism/LineSearch.cpp base/globalization_mechanism/TrustLineSearch.cpp \
	base/local_approximation/QP.cpp base/local_approximation/LP.cpp \
	base/local_approximation/LocalApproximation.cpp base/local_approximation/QPApproximation.cpp \
	base/local_approximation/LocalSolution.cpp base/Utils.cpp \
	solvers/BQPDSolver.cpp
OBJ_FILES = $(subst .cpp,.o,$(SRC_FILES)) solvers/wdotd.o
EXEC_NAME = argonot

# base/globalization_strategy/Tube.cpp base/globalization_strategy/TubeGlobalizationStrategy.cpp
# base/direction_computation/LPDirectionComputation.cpp \

DEFINES =# -D BQPD # https://stackoverflow.com/questions/10337248/what-does-the-g-d-flag-do
LIBS = $(AMPL_DIR)/amplsolver.a $(BQPD_DIR)/libbqpd.a

.cpp.o:
	$(CPP) $(CFLAGS) $(DEFINES) $(INCLUDES) -c $< -o $@ ;

.c.o:
	$(CC) -c $(FFLAGS) $(CFLAGS) $(DEFINES) $(INCLUDES) $*.c -o $@

.f.o:
	$(F77) -c $(FFLAGS) $< -o $@

.f90.o:
	$(F77) -I $(FFLAGS) -c $< -o $@ ;

all: $(EXEC_NAME) $(AMPL_DIR)/amplsolver.a

$(EXEC_NAME): $(OBJ_FILES) $(LIBS)
	$(CPP) $(LFLAGS) -o $(EXEC_NAME) $(OBJ_FILES) $(LIBS) $(LIBLINALG) $(L) $(LIBCOMPIL)

$(AMPL_DIR)/amplsolver.a:
	cd $(AMPL_DIR); make amplsolver.a

clean: cleanhere

cleanhere:
	rm -f $(EXEC_NAME) $(OBJ_FILES) *~

xsum.out: $SRC_FILES
	xsum $SRC_FILES >xsum1.out
	cmp xsum0.out xsum1.out && mv xsum1.out xsum.out || diff xsum[01].out
