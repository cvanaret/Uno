#  NLP written by GAMS Convert at 06/20/02 11:29:54
#  
#  Equation counts
#     Total       E       G       L       N       X
#         3       2       1       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         5       5       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        12       6       6       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 := 50, >= 50, <= 200;
var x2 := 37.5, >= 37.5, <= 150;
var x3 := 45, >= 45, <= 180;
var x4;

minimize obj: 0.00533*x1^2 + 11.669*x1 + 0.00889*x2^2 + 10.333*x2 + 0.00741*x3^
              2 + 10.833*x3 + 653.1;

subject to

e2:  - (0.01*(0.0676*x1*x1 + 0.00953*x1*x2 - 0.00507*x1*x3 + 0.00953*x2*x1 + 
    0.0521*x2*x2 + 0.00901*x2*x3 - 0.00507*x3*x1 + 0.00901*x3*x2 + 0.0294*x3*x3
    ) - 0.000766*x1 - 3.42e-5*x2 + 0.000189*x3) + x4 = 0.040357;

e3:    x1 + x2 + x3 - x4 >= 210;
