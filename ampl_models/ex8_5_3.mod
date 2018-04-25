#  NLP written by GAMS Convert at 07/05/02 23:00:40
#  
#  Equation counts
#     Total       E       G       L       N       X       C
#         7       5       2       0       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         6       6       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        20      10      10       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x2 := 0.5, >= 0;
var x3 := 0.5, >= 0;
var x4 := 2, >= 0;
var x5;
var x6 := 1, >= 0;

minimize obj: x2*log(x2) + x3*log(x3) - log(x4 - x6) + x4 - x5*log(1 + x6/x4)/
              x6 + 5.0464317551216*x2 + 0.366877055769689*x3 - 1;

subject to

e2: x4^3 - x4^2 + (-x6^2 - x6 + x5)*x4 - x5*x6 = 0;

e3:  - (1.04633*x2*x2 + 0.579822*x2*x3 + 0.579822*x3*x2 + 0.379615*x3*x3) + x5
     = 0;

e4:  - 0.0771517*x2 - 0.0765784*x3 + x6 = 0;

e5:    x2 + x3 = 1;

e6:    x4 - x6 >= 0;

e7:    x2 >= 0.05;
