#  NLP written by GAMS Convert at 06/20/02 11:27:32
#  
#  Equation counts
#     Total       E       G       L       N       X
#         4       2       2       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         5       5       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        17      13       4       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x2 := 0.685244910300343, >= 0;
var x3 := 0.0126990526103601, >= 0;
var x4 := 0.302056037089293, >= 0;
var x5 >= 0;

minimize obj:    24.55*x2 + 26.75*x3 + 39*x4 + 40.5*x5;

subject to

e2:    x2 + x3 + x4 + x5 = 1;

e3: 12*x2 - 1.645*sqrt(0.28*x2^2 + 0.19*x3^2 + 20.5*x4^2 + 0.62*x5^2) + 11.9*x3
     + 41.8*x4 + 52.1*x5 >= 21;

e4:    2.3*x2 + 5.6*x3 + 11.1*x4 + 1.3*x5 >= 5;
