#  NLP written by GAMS Convert at 06/20/02 11:42:10
#  
#  Equation counts
#     Total       E       G       L       N       X
#        11      11       0       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        25      25       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        73      49      24       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 := 6, >= 0, <= 100;
var x2 := 2, >= 0, <= 100;
var x3 >= 0, <= 100;
var x4 >= 0, <= 100;
var x5 >= 0, <= 100;
var x6 := 3, >= 0, <= 100;
var x7 >= 0, <= 100;
var x8 := 21, >= 0, <= 100;
var x9 := 20, >= 0, <= 100;
var x10 >= 0, <= 100;
var x11 >= 0, <= 100;
var x12 >= 0, <= 100;
var x13 := 24, >= 0, <= 100;
var x14 >= 0, <= 100;
var x15 >= 0, <= 100;
var x16 >= 0, <= 100;
var x17 := 3, >= 0, <= 100;
var x18 >= 0, <= 100;
var x19 := 13, >= 0, <= 100;
var x20 >= 0, <= 100;
var x21 >= 0, <= 100;
var x22 := 12, >= 0, <= 100;
var x23 >= 0, <= 100;
var x24 >= 0, <= 100;

minimize obj: 300*x1 - 7*x1*x1 - 4*x2*x2 + 270*x2 - 6*x3*x3 + 460*x3 - 8*x4*x4
               + 800*x4 - 12*x5*x5 + 740*x5 - 9*x6*x6 + 600*x6 - 14*x7*x7 + 540
              *x7 - 7*x8*x8 + 380*x8 - 13*x9*x9 + 300*x9 - 12*x10*x10 + 490*x10
               - 8*x11*x11 + 380*x11 - 4*x12*x12 + 760*x12 - 7*x13*x13 + 430*
              x13 - 9*x14*x14 + 250*x14 - 16*x15*x15 + 390*x15 - 8*x16*x16 + 
              600*x16 - 4*x17*x17 + 210*x17 - 10*x18*x18 + 830*x18 - 21*x19*x19
               + 470*x19 - 13*x20*x20 + 680*x20 - 17*x21*x21 + 360*x21 - 9*x22*
              x22 + 290*x22 - 8*x23*x23 + 400*x23 - 4*x24*x24 + 310*x24;

subject to

e2:    x1 + x2 + x3 + x4 = 8;

e3:    x5 + x6 + x7 + x8 = 24;

e4:    x9 + x10 + x11 + x12 = 20;

e5:    x13 + x14 + x15 + x16 = 24;

e6:    x17 + x18 + x19 + x20 = 16;

e7:    x21 + x22 + x23 + x24 = 12;

e8:    x1 + x5 + x9 + x13 + x17 + x21 = 29;

e9:    x2 + x6 + x10 + x14 + x18 + x22 = 41;

e10:    x3 + x7 + x11 + x15 + x19 + x23 = 13;

e11:    x4 + x8 + x12 + x16 + x20 + x24 = 21;
