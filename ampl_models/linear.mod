#  DNLP written by GAMS Convert at 06/20/02 12:27:16
#  
#  Equation counts
#     Total       E       G       L       N       X
#        21      21       0       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        25      25       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#       121     101      20       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x2 >= -100, <= 100;
var x3 >= -100, <= 100;
var x4 := -92, >= -100, <= 100;
var x5 := -94, >= -100, <= 100;
var x6 >= -100, <= 100;
var x7 := -94, >= -100, <= 100;
var x8 := -96, >= -100, <= 100;
var x9 := -83, >= -100, <= 100;
var x10 := -90, >= -100, <= 100;
var x11 := -93, >= -100, <= 100;
var x12 >= -100, <= 100;
var x13 >= -100, <= 100;
var x14 >= -100, <= 100;
var x15 >= -100, <= 100;
var x16 >= -100, <= 100;
var x17 >= -100, <= 100;
var x18 := -84, >= -100, <= 100;
var x19 := -83, >= -100, <= 100;
var x20 := -92, >= -100, <= 100;
var x21 >= -100, <= 100;
var x22 := 1;
var x23 := 1;
var x24 := 1;
var x25 := 1;

minimize obj: abs(x2) + abs(x3) + abs(x4) + abs(x5) + abs(x6) + abs(x7) + abs(
              x8) + abs(x9) + abs(x10) + abs(x11) + abs(x12) + abs(x13) + abs(
              x14) + abs(x15) + abs(x16) + abs(x17) + abs(x18) + abs(x19) + 
              abs(x20) + abs(x21);

subject to

e1:    x2 + x22 + 85*x23 + 76*x24 + 44*x25 = 99;

e2:    x3 + x22 + 82*x23 + 78*x24 + 42*x25 = 93;

e3:    x4 + x22 + 75*x23 + 73*x24 + 42*x25 = 99;

e4:    x5 + x22 + 74*x23 + 72*x24 + 44*x25 = 97;

e5:    x6 + x22 + 76*x23 + 73*x24 + 43*x25 = 90;

e6:    x7 + x22 + 74*x23 + 69*x24 + 46*x25 = 96;

e7:    x8 + x22 + 73*x23 + 69*x24 + 46*x25 = 93;

e8:    x9 + x22 + 96*x23 + 80*x24 + 36*x25 = 130;

e9:    x10 + x22 + 93*x23 + 78*x24 + 36*x25 = 118;

e10:    x11 + x22 + 70*x23 + 73*x24 + 37*x25 = 88;

e11:    x12 + x22 + 82*x23 + 71*x24 + 46*x25 = 89;

e12:    x13 + x22 + 80*x23 + 72*x24 + 45*x25 = 93;

e13:    x14 + x22 + 77*x23 + 76*x24 + 42*x25 = 94;

e14:    x15 + x22 + 67*x23 + 76*x24 + 50*x25 = 75;

e15:    x16 + x22 + 82*x23 + 70*x24 + 48*x25 = 84;

e16:    x17 + x22 + 76*x23 + 76*x24 + 41*x25 = 91;

e17:    x18 + x22 + 74*x23 + 78*x24 + 31*x25 = 100;

e18:    x19 + x22 + 71*x23 + 80*x24 + 29*x25 = 98;

e19:    x20 + x22 + 70*x23 + 83*x24 + 39*x25 = 101;

e20:    x21 + x22 + 64*x23 + 79*x24 + 38*x25 = 80;
