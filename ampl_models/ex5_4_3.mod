#  NLP written by GAMS Convert at 06/20/02 11:50:34
#  
#  Equation counts
#     Total       E       G       L       N       X
#        14      14       0       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        17      17       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        43      25      18       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 := 10, >= 10, <= 350;
var x2 := 10, >= 10, <= 350;
var x3 := 10, >= 10, <= 200;
var x4 := 10, >= 10, <= 200;
var x5 >= 0, <= 10;
var x6 >= 0, <= 10;
var x7 >= 0, <= 10;
var x8 >= 0, <= 10;
var x9 >= 0, <= 10;
var x10 >= 0, <= 10;
var x11 >= 0, <= 10;
var x12 >= 0, <= 10;
var x13 := 150, >= 150, <= 310;
var x14 := 150, >= 150, <= 310;
var x15 := 150, >= 150, <= 310;
var x16 := 150, >= 150, <= 310;

minimize obj: 1300*(1000/(0.0333333333333333*x1*x2 + 0.166666666666667*x1 + 
              0.166666666666667*x2))^0.6 + 1300*(600/(0.0333333333333333*x3*x4
               + 0.166666666666667*x3 + 0.166666666666667*x4))^0.6;

subject to

e1:    x5 + x9 = 10;

e2:    x5 - x6 + x11 = 0;

e3:    x7 + x9 - x10 = 0;

e4:  - x6 + x7 + x8 = 0;

e5:  - x10 + x11 + x12 = 0;

e6: x16*x11 - x13*x6 + 150*x5 = 0;

e7: x15*x7 - x14*x10 + 150*x9 = 0;

e8: x6*x15 - x6*x13 = 1000;

e9: x10*x16 - x10*x14 = 600;

e10:    x1 + x15 = 500;

e11:    x2 + x13 = 250;

e12:    x3 + x16 = 350;

e13:    x4 + x14 = 200;
