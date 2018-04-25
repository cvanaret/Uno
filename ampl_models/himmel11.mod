#  NLP written by GAMS Convert at 06/20/02 12:23:00
#  
#  Equation counts
#     Total       E       G       L       N       X
#         5       4       1       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        10      10       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        27      11      16       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 >= 0, <= 92;
var x2 := 90, >= 90, <= 110;
var x3 := 20, >= 20, <= 25;
var x4 >= 0;
var x5 := 78.62, >= 78, <= 102;
var x6 := 33.44, >= 33, <= 45;
var x7 := 31.07, >= 27, <= 45;
var x8 := 44.18, >= 27, <= 45;
var x9 := 35.22, >= 27, <= 45;

minimize obj: 5.3578547*x7^2 + 0.8356891*x5*x9 + 37.293239*x5 + 5000*x4
               - 40792.141;

subject to

e1:    5*x4 - x5 + 7*x7 - x9 >= 0;

e2:  - (0.0056858*x6*x9 + 0.0006262*x5*x8 - 0.0022053*x7*x9) + x1 + 2*x4
     = 85.334407;

e3:  - (0.0071317*x6*x9 + 0.0029955*x5*x6 + 0.0021813*x7^2) + x2 = 80.51249;

e4:  - (0.0047026*x7*x9 + 0.0012547*x5*x7 + 0.0019085*x7*x8) + x3 + 4*x4
     = 9.300961;
