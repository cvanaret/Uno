#  NLP written by GAMS Convert at 06/20/02 11:33:47
#  
#  Equation counts
#     Total       E       G       L       N       X
#         5       1       0       4       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         4       4       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        14       6       8       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 := 5.49E-6, >= 5.49E-6, <= 4.553;
var x2 := 0.0021961, >= 0.0021961, <= 18.21;
var x3;

minimize obj:    x3;

subject to

e2: 10000*x1*x2 - x3 <= 1;

e3:  - 10000*x1*x2 - x3 <= -1;

e4: exp(-x1) + exp(-x2) - x3 <= 1.001;

e5: (-exp(-x1)) - exp(-x2) - x3 <= -1.001;
