#  NLP written by GAMS Convert at 06/20/02 11:55:50
#  
#  Equation counts
#     Total       E       G       L       N       X
#         6       5       0       1       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         7       7       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        17       7      10       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 >= 0, <= 1;
var x2 >= 0, <= 1;
var x3 >= 0, <= 1;
var x4 >= 0, <= 1;
var x5 := 1E-5, >= 1E-5, <= 16;
var x6 := 1E-5, >= 1E-5, <= 16;

minimize obj:  - x4;

subject to

e2: 0.09755988*x1*x5 + x1 = 1;

e3: 0.0965842812*x2*x6 + x2 - x1 = 0;

e4: 0.0391908*x3*x5 + x3 + x1 = 1;

e5: 0.03527172*x4*x6 + x4 - x1 + x2 - x3 = 0;

e6: x5**0.5 + x6**0.5 <= 4;
