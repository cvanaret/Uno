#  NLP written by GAMS Convert at 06/20/02 11:21:07
#  
#  Equation counts
#     Total       E       G       L       N       X
#        13      10       1       2       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        14      14       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        41      13      28       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 := 6, >= 1, <= 16;
var x2 := 5, >= 1, <= 16;
var x3 := 6, >= 1, <= 16;
var x4 := 3, >= 1, <= 16;
var x6 := 1000, >= 1, <= 1000;
var x7 := 1.6, >= 0.0001;
var x8 := 0.3, >= 0.0001;
var x9 := 1, >= 1;
var x10 := 50, <= 50;
var x11 := 600, >= 100;
var x12 := 1, >= 1;
var x13 := 0.0001, >= 0.0001;
var x14 := 0.01, >= 0.01;

minimize obj:    x7 + x8;

subject to

e2:  - 1.42857142857143*x4*x6 + 10000*x8 = 0;

e3: 10*x7*x9 - 0.00968946189201592*x3*(x1**4 - x2**4) = 0;

e4: 143.3076*x10*x4 - 10000*x7 = 0;

e5: 3.1415927*x6*(0.001*x9)**3 - 6e-6*x3*x4*x13 = 0;

e6: 101000*x12*x13 - 1.57079635*x6*x14 = 0;

e7: log10(0.8 + 8.112*x3) - 10964781961.4318*x11**(-3.55) = 0;

e8:  - 0.5*x10 + x11 = 560;

e9:    x1 - x2 >= 0;

e10: 0.0307*x4^2 - 0.3864*(0.0062831854*x1*x9)^2*x6 <= 0;

e11:    101000*x12 - 15707.9635*x14 <= 0;

e12:  - (log(x1) - log(x2)) + x13 = 0;

e13:  - (x1^2 - x2^2) + x14 = 0;
