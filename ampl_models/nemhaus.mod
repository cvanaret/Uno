#  NLP written by GAMS Convert at 07/09/02 21:33:15
#  
#  Equation counts
#     Total       E       G       L       N       X       C
#         6       6       0       0       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         6       6       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        11       6       5       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x2 >= 0;
var x3 >= 0;
var x4 >= 0;
var x5 >= 0;
var x6 >= 0;

minimize obj: 2*x2*x4 + 4*x2*x5 + 3*x2*x6 + 6*x3*x4 + 2*x3*x5 + 3*x3*x6 + 5*x4*
              x5 + 3*x4*x6 + 3*x5*x6;

subject to

e2:    x2 = 1;

e3:    x3 = 1;

e4:    x4 = 1;

e5:    x5 = 1;

e6:    x6 = 1;
