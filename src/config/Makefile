.SUFFIXES: .c .f .o .f90

# specify local compilers and compiler flags
CC     = gcc
CPP    = g++
F77    = gfortran
FFLAGS = -g 
CFLAGS = -Wall -pedantic -Wno-write-strings -std=c++11 -fpermissive
LFLAGS = 

INCLUDES =

SRC_FILES = main.cpp ConfigurationFile.cpp
OBJ_FILES = $(subst .cpp,.o,$(SRC_FILES))
EXEC_NAME = config

DEFINES =
LIBS =

.cpp.o:
	$(CPP) $(CFLAGS) $(DEFINES) $(INCLUDES) -c $< -o $@ ;

all: config

config: $(OBJ_FILES)
	$(CPP) $(LFLAGS) -o $(EXEC_NAME) $(OBJ_FILES) $(LIBS)

clean:
	rm -f $(EXEC_NAME) $(OBJ_FILES) *~

xsum.out: $SRC_FILES
	xsum $SRC_FILES >xsum1.out
	cmp xsum0.out xsum1.out && mv xsum1.out xsum.out || diff xsum[01].out
