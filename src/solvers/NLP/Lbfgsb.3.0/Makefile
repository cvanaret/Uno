FC = gfortran

FFLAGS = -O -Wall -fbounds-check -g -Wno-uninitialized 

DRIVER1_77 = driver1.f
DRIVER2_77 = driver2.f
DRIVER3_77 = driver3.f

DRIVER1_90 = driver1.f90
DRIVER2_90 = driver2.f90
DRIVER3_90 = driver3.f90
DRIVERSL1_90 = driverSL1.f90
DRIVERCR_90 = driverChainRosen.f90

LBFGSB  = lbfgsb.f
LINPACK = linpack.f
BLAS    = blas.f
TIMER   = timer.f

all :  lbfgsb_77_1 lbfgsb_77_2 lbfgsb_77_3 lbfgsb_90_1 lbfgsb_90_2 lbfgsb_90_3 lbfgsb_90_SL1 lbfgsb_90_CR  


lbfgsb_77_1 : $(DRIVER1_77) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(FC) $(FFLAGS) $(DRIVER1_77) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) -o x.lbfgsb_77_1

lbfgsb_77_2 : $(DRIVER2_77) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(FC) $(FFLAGS) $(DRIVER2_77) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) -o x.lbfgsb_77_2

lbfgsb_77_3 : $(DRIVER3_77) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(FC) $(FFLAGS) $(DRIVER3_77) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) -o x.lbfgsb_77_3

lbfgsb_90_1 : $(DRIVER1_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(FC) $(FFLAGS) $(DRIVER1_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) -o x.lbfgsb_90_1

lbfgsb_90_2 : $(DRIVER2_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(FC) $(FFLAGS) $(DRIVER2_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) -o x.lbfgsb_90_2

lbfgsb_90_3 : $(DRIVER3_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(FC) $(FFLAGS) $(DRIVER3_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) -o x.lbfgsb_90_3

lbfgsb_90_SL1 : $(DRIVERSL1_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(FC) $(FFLAGS) $(DRIVERSL1_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) -o x.lbfgsb_90_SL1

lbfgsb_90_CR : $(DRIVERCR_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(FC) $(FFLAGS) $(DRIVERCR_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) -o x.lbfgsb_90_CR

clean:
	rm -f *.o x.lbfgsb* *~
