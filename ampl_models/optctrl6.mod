#
#**************************
# SET UP THE INITIAL DATA *
#**************************
#   Problem :
#   *********
#   An optimal control problem
#   This problem is a modification of OPTCNTRL.  All bound constraints
#   were removed.  At the solution, the bounds on U1-U9 are active,
#   so a penalty term was added to the objective:
#
#        mu * ||u(i) - active_bound||
#
#   The value of mu (TDP_mu) was chosen to be very large, making the
#   Hessian of the objective ill-conditioned.
#   Source:
#   B. Murtagh and M. Saunders,
#   Mathematical Programming studies 16, pp 84-117,
#   (example 5.11)
#   SIF input: Nick Gould, June 1990.
#              modified by T. Plantagena, December 1992.
#   classification QQR2-AN-V-V
#   useful parameters
#IE T                   10
#   penalty parameter mu
#   Least square problems are bounded below by zero
#   Solution
param t := 40;
param tm1 := -1 + (40);
param tdp_mu := 1000000.0;
param tp1 := 1 + (39);

var x0 >= 10.0 ,  <= 10.0;
var y0 >= 0.0 ,  <= 0.0;
var x1;
var y1 := -1.0;
var x2;
var y2 := -1.0;
var x3;
var y3 := -1.0;
var x4;
var y4 := -1.0;
var x5;
var y5 := -1.0;
var x6;
var y6 := -1.0;
var x7;
var y7 := -1.0;
var x8;
var y8 := -1.0;
var x9;
var y9 := -1.0;
var x10;
var y10 := -1.0;
var x11;
var y11 := -1.0;
var x12;
var y12 := -1.0;
var x13;
var y13 := -1.0;
var x14;
var y14 := -1.0;
var x15;
var y15 := -1.0;
var x16;
var y16 := -1.0;
var x17;
var y17 := -1.0;
var x18;
var y18 := -1.0;
var x19;
var y19 := -1.0;
var x20;
var y20 := -1.0;
var x21;
var y21 := -1.0;
var x22;
var y22 := -1.0;
var x23;
var y23 := -1.0;
var x24;
var y24 := -1.0;
var x25;
var y25 := -1.0;
var x26;
var y26 := -1.0;
var x27;
var y27 := -1.0;
var x28;
var y28 := -1.0;
var x29;
var y29 := -1.0;
var x30;
var y30 := -1.0;
var x31;
var y31 := -1.0;
var x32;
var y32 := -1.0;
var x33;
var y33 := -1.0;
var x34;
var y34 := -1.0;
var x35;
var y35 := -1.0;
var x36;
var y36 := -1.0;
var x37;
var y37 := -1.0;
var x38;
var y38 := -1.0;
var x39;
var y39 := -1.0;
var x40;
var y40 >= 0.0 ,  <= 0.0;
var u0;
var u1;
var u2;
var u3;
var u4;
var u5;
var u6;
var u7;
var u8;
var u9;
var u10;
var u11;
var u12;
var u13;
var u14;
var u15;
var u16;
var u17;
var u18;
var u19;
var u20;
var u21;
var u22;
var u23;
var u24;
var u25;
var u26;
var u27;
var u28;
var u29;
var u30;
var u31;
var u32;
var u33;
var u34;
var u35;
var u36;
var u37;
var u38;
var u39;

minimize obj:
	0.5*x0 * x0 + 0.5*x1 * x1 + 0.5*x2 * x2 + 0.5*x3 * x3 + 0.5*x4 * x4 + 0.5*x5 * x5 + 0.5*x6 * x6 + 0.5*x7 * x7 + 0.5*x8 * x8 + 0.5*x9 * x9 + 0.5*x10 * x10 + 0.5*x11 * x11 + 0.5*x12 * x12 + 0.5*x13 * x13 + 0.5*x14 * x14 + 0.5*x15 * x15 + 0.5*x16 * x16 + 0.5*x17 * x17 + 0.5*x18 * x18 + 0.5*x19 * x19 + 0.5*x20 * x20 + 0.5*x21 * x21 + 0.5*x22 * x22 + 0.5*x23 * x23 + 0.5*x24 * x24 + 0.5*x25 * x25 + 0.5*x26 * x26 + 0.5*x27 * x27 + 0.5*x28 * x28 + 0.5*x29 * x29 + 0.5*x30 * x30 + 0.5*x31 * x31 + 0.5*x32 * x32 + 0.5*x33 * x33 + 0.5*x34 * x34 + 0.5*x35 * x35 + 0.5*x36 * x36 + 0.5*x37 * x37 + 0.5*x38 * x38 + 0.5*x39 * x39 + 0.5*x40 * x40 + 1000000.0*(u1-0.2) * (u1-0.2) + 1000000.0*(u2-0.2) * (u2-0.2) + 1000000.0*(u3-0.2) * (u3-0.2) + 1000000.0*(u4-0.2) * (u4-0.2) + 1000000.0*(u5-0.2) * (u5-0.2) + 1000000.0*(u6-0.2) * (u6-0.2) + 1000000.0*(u7-0.2) * (u7-0.2) + 1000000.0*(u8-0.2) * (u8-0.2) + 1000000.0*(u9-0.2) * (u9-0.2) + 1000000.0*(u10-0.2) * (u10-0.2) + 1000000.0*(u11-0.2) * (u11-0.2) + 1000000.0*(u12-0.2) * (u12-0.2) + 1000000.0*(u13-0.2) * (u13-0.2) + 1000000.0*(u14-0.2) * (u14-0.2) + 1000000.0*(u15-0.2) * (u15-0.2) + 1000000.0*(u16-0.2) * (u16-0.2) + 1000000.0*(u17-0.2) * (u17-0.2) + 1000000.0*(u18-0.2) * (u18-0.2) + 1000000.0*(u19-0.2) * (u19-0.2) + 1000000.0*(u20-0.2) * (u20-0.2) + 1000000.0*(u21-0.2) * (u21-0.2) + 1000000.0*(u22-0.2) * (u22-0.2) + 1000000.0*(u23-0.2) * (u23-0.2) + 1000000.0*(u24-0.2) * (u24-0.2) + 1000000.0*(u25-0.2) * (u25-0.2) + 1000000.0*(u26-0.2) * (u26-0.2) + 1000000.0*(u27-0.2) * (u27-0.2) + 1000000.0*(u28-0.2) * (u28-0.2) + 1000000.0*(u29-0.2) * (u29-0.2) + 1000000.0*(u30-0.2) * (u30-0.2) + 1000000.0*(u31-0.2) * (u31-0.2) + 1000000.0*(u32-0.2) * (u32-0.2) + 1000000.0*(u33-0.2) * (u33-0.2) + 1000000.0*(u34-0.2) * (u34-0.2) + 1000000.0*(u35-0.2) * (u35-0.2) + 1000000.0*(u36-0.2) * (u36-0.2) + 1000000.0*(u37-0.2) * (u37-0.2) + 1000000.0*(u38-0.2) * (u38-0.2) + 1000000.0*(u39-0.2) * (u39-0.2);

subject to obj_bnd:
	0.0 <= 0.5*x0 * x0 + 0.5*x1 * x1 + 0.5*x2 * x2 + 0.5*x3 * x3 + 0.5*x4 * x4 + 0.5*x5 * x5 + 0.5*x6 * x6 + 0.5*x7 * x7 + 0.5*x8 * x8 + 0.5*x9 * x9 + 0.5*x10 * x10 + 0.5*x11 * x11 + 0.5*x12 * x12 + 0.5*x13 * x13 + 0.5*x14 * x14 + 0.5*x15 * x15 + 0.5*x16 * x16 + 0.5*x17 * x17 + 0.5*x18 * x18 + 0.5*x19 * x19 + 0.5*x20 * x20 + 0.5*x21 * x21 + 0.5*x22 * x22 + 0.5*x23 * x23 + 0.5*x24 * x24 + 0.5*x25 * x25 + 0.5*x26 * x26 + 0.5*x27 * x27 + 0.5*x28 * x28 + 0.5*x29 * x29 + 0.5*x30 * x30 + 0.5*x31 * x31 + 0.5*x32 * x32 + 0.5*x33 * x33 + 0.5*x34 * x34 + 0.5*x35 * x35 + 0.5*x36 * x36 + 0.5*x37 * x37 + 0.5*x38 * x38 + 0.5*x39 * x39 + 0.5*x40 * x40 + 1000000.0*(u1-0.2) * (u1-0.2) + 1000000.0*(u2-0.2) * (u2-0.2) + 1000000.0*(u3-0.2) * (u3-0.2) + 1000000.0*(u4-0.2) * (u4-0.2) + 1000000.0*(u5-0.2) * (u5-0.2) + 1000000.0*(u6-0.2) * (u6-0.2) + 1000000.0*(u7-0.2) * (u7-0.2) + 1000000.0*(u8-0.2) * (u8-0.2) + 1000000.0*(u9-0.2) * (u9-0.2) + 1000000.0*(u10-0.2) * (u10-0.2) + 1000000.0*(u11-0.2) * (u11-0.2) + 1000000.0*(u12-0.2) * (u12-0.2) + 1000000.0*(u13-0.2) * (u13-0.2) + 1000000.0*(u14-0.2) * (u14-0.2) + 1000000.0*(u15-0.2) * (u15-0.2) + 1000000.0*(u16-0.2) * (u16-0.2) + 1000000.0*(u17-0.2) * (u17-0.2) + 1000000.0*(u18-0.2) * (u18-0.2) + 1000000.0*(u19-0.2) * (u19-0.2) + 1000000.0*(u20-0.2) * (u20-0.2) + 1000000.0*(u21-0.2) * (u21-0.2) + 1000000.0*(u22-0.2) * (u22-0.2) + 1000000.0*(u23-0.2) * (u23-0.2) + 1000000.0*(u24-0.2) * (u24-0.2) + 1000000.0*(u25-0.2) * (u25-0.2) + 1000000.0*(u26-0.2) * (u26-0.2) + 1000000.0*(u27-0.2) * (u27-0.2) + 1000000.0*(u28-0.2) * (u28-0.2) + 1000000.0*(u29-0.2) * (u29-0.2) + 1000000.0*(u30-0.2) * (u30-0.2) + 1000000.0*(u31-0.2) * (u31-0.2) + 1000000.0*(u32-0.2) * (u32-0.2) + 1000000.0*(u33-0.2) * (u33-0.2) + 1000000.0*(u34-0.2) * (u34-0.2) + 1000000.0*(u35-0.2) * (u35-0.2) + 1000000.0*(u36-0.2) * (u36-0.2) + 1000000.0*(u37-0.2) * (u37-0.2) + 1000000.0*(u38-0.2) * (u38-0.2) + 1000000.0*(u39-0.2) * (u39-0.2);
subject to b0:
	x1 - x0 - 0.2*y0 = 0;
subject to c0:
	0.01*y0 * y0 + y1 - y0 + 0.0040*x0 - 0.2*u0 = 0;
subject to b1:
	x2 - x1 - 0.2*y1 = 0;
subject to c1:
	0.01*y1 * y1 + y2 - y1 + 0.0040*x1 - 0.2*u1 = 0;
subject to b2:
	x3 - x2 - 0.2*y2 = 0;
subject to c2:
	0.01*y2 * y2 + y3 - y2 + 0.0040*x2 - 0.2*u2 = 0;
subject to b3:
	x4 - x3 - 0.2*y3 = 0;
subject to c3:
	0.01*y3 * y3 + y4 - y3 + 0.0040*x3 - 0.2*u3 = 0;
subject to b4:
	x5 - x4 - 0.2*y4 = 0;
subject to c4:
	0.01*y4 * y4 + y5 - y4 + 0.0040*x4 - 0.2*u4 = 0;
subject to b5:
	x6 - x5 - 0.2*y5 = 0;
subject to c5:
	0.01*y5 * y5 + y6 - y5 + 0.0040*x5 - 0.2*u5 = 0;
subject to b6:
	x7 - x6 - 0.2*y6 = 0;
subject to c6:
	0.01*y6 * y6 + y7 - y6 + 0.0040*x6 - 0.2*u6 = 0;
subject to b7:
	x8 - x7 - 0.2*y7 = 0;
subject to c7:
	0.01*y7 * y7 + y8 - y7 + 0.0040*x7 - 0.2*u7 = 0;
subject to b8:
	x9 - x8 - 0.2*y8 = 0;
subject to c8:
	0.01*y8 * y8 + y9 - y8 + 0.0040*x8 - 0.2*u8 = 0;
subject to b9:
	x10 - x9 - 0.2*y9 = 0;
subject to c9:
	0.01*y9 * y9 + y10 - y9 + 0.0040*x9 - 0.2*u9 = 0;
subject to b10:
	x11 - x10 - 0.2*y10 = 0;
subject to c10:
	0.01*y10 * y10 + y11 - y10 + 0.0040*x10 - 0.2*u10 = 0;
subject to b11:
	x12 - x11 - 0.2*y11 = 0;
subject to c11:
	0.01*y11 * y11 + y12 - y11 + 0.0040*x11 - 0.2*u11 = 0;
subject to b12:
	x13 - x12 - 0.2*y12 = 0;
subject to c12:
	0.01*y12 * y12 + y13 - y12 + 0.0040*x12 - 0.2*u12 = 0;
subject to b13:
	x14 - x13 - 0.2*y13 = 0;
subject to c13:
	0.01*y13 * y13 + y14 - y13 + 0.0040*x13 - 0.2*u13 = 0;
subject to b14:
	x15 - x14 - 0.2*y14 = 0;
subject to c14:
	0.01*y14 * y14 + y15 - y14 + 0.0040*x14 - 0.2*u14 = 0;
subject to b15:
	x16 - x15 - 0.2*y15 = 0;
subject to c15:
	0.01*y15 * y15 + y16 - y15 + 0.0040*x15 - 0.2*u15 = 0;
subject to b16:
	x17 - x16 - 0.2*y16 = 0;
subject to c16:
	0.01*y16 * y16 + y17 - y16 + 0.0040*x16 - 0.2*u16 = 0;
subject to b17:
	x18 - x17 - 0.2*y17 = 0;
subject to c17:
	0.01*y17 * y17 + y18 - y17 + 0.0040*x17 - 0.2*u17 = 0;
subject to b18:
	x19 - x18 - 0.2*y18 = 0;
subject to c18:
	0.01*y18 * y18 + y19 - y18 + 0.0040*x18 - 0.2*u18 = 0;
subject to b19:
	x20 - x19 - 0.2*y19 = 0;
subject to c19:
	0.01*y19 * y19 + y20 - y19 + 0.0040*x19 - 0.2*u19 = 0;
subject to b20:
	x21 - x20 - 0.2*y20 = 0;
subject to c20:
	0.01*y20 * y20 + y21 - y20 + 0.0040*x20 - 0.2*u20 = 0;
subject to b21:
	x22 - x21 - 0.2*y21 = 0;
subject to c21:
	0.01*y21 * y21 + y22 - y21 + 0.0040*x21 - 0.2*u21 = 0;
subject to b22:
	x23 - x22 - 0.2*y22 = 0;
subject to c22:
	0.01*y22 * y22 + y23 - y22 + 0.0040*x22 - 0.2*u22 = 0;
subject to b23:
	x24 - x23 - 0.2*y23 = 0;
subject to c23:
	0.01*y23 * y23 + y24 - y23 + 0.0040*x23 - 0.2*u23 = 0;
subject to b24:
	x25 - x24 - 0.2*y24 = 0;
subject to c24:
	0.01*y24 * y24 + y25 - y24 + 0.0040*x24 - 0.2*u24 = 0;
subject to b25:
	x26 - x25 - 0.2*y25 = 0;
subject to c25:
	0.01*y25 * y25 + y26 - y25 + 0.0040*x25 - 0.2*u25 = 0;
subject to b26:
	x27 - x26 - 0.2*y26 = 0;
subject to c26:
	0.01*y26 * y26 + y27 - y26 + 0.0040*x26 - 0.2*u26 = 0;
subject to b27:
	x28 - x27 - 0.2*y27 = 0;
subject to c27:
	0.01*y27 * y27 + y28 - y27 + 0.0040*x27 - 0.2*u27 = 0;
subject to b28:
	x29 - x28 - 0.2*y28 = 0;
subject to c28:
	0.01*y28 * y28 + y29 - y28 + 0.0040*x28 - 0.2*u28 = 0;
subject to b29:
	x30 - x29 - 0.2*y29 = 0;
subject to c29:
	0.01*y29 * y29 + y30 - y29 + 0.0040*x29 - 0.2*u29 = 0;
subject to b30:
	x31 - x30 - 0.2*y30 = 0;
subject to c30:
	0.01*y30 * y30 + y31 - y30 + 0.0040*x30 - 0.2*u30 = 0;
subject to b31:
	x32 - x31 - 0.2*y31 = 0;
subject to c31:
	0.01*y31 * y31 + y32 - y31 + 0.0040*x31 - 0.2*u31 = 0;
subject to b32:
	x33 - x32 - 0.2*y32 = 0;
subject to c32:
	0.01*y32 * y32 + y33 - y32 + 0.0040*x32 - 0.2*u32 = 0;
subject to b33:
	x34 - x33 - 0.2*y33 = 0;
subject to c33:
	0.01*y33 * y33 + y34 - y33 + 0.0040*x33 - 0.2*u33 = 0;
subject to b34:
	x35 - x34 - 0.2*y34 = 0;
subject to c34:
	0.01*y34 * y34 + y35 - y34 + 0.0040*x34 - 0.2*u34 = 0;
subject to b35:
	x36 - x35 - 0.2*y35 = 0;
subject to c35:
	0.01*y35 * y35 + y36 - y35 + 0.0040*x35 - 0.2*u35 = 0;
subject to b36:
	x37 - x36 - 0.2*y36 = 0;
subject to c36:
	0.01*y36 * y36 + y37 - y36 + 0.0040*x36 - 0.2*u36 = 0;
subject to b37:
	x38 - x37 - 0.2*y37 = 0;
subject to c37:
	0.01*y37 * y37 + y38 - y37 + 0.0040*x37 - 0.2*u37 = 0;
subject to b38:
	x39 - x38 - 0.2*y38 = 0;
subject to c38:
	0.01*y38 * y38 + y39 - y38 + 0.0040*x38 - 0.2*u38 = 0;
subject to b39:
	x40 - x39 - 0.2*y39 = 0;
subject to c39:
	0.01*y39 * y39 + y40 - y39 + 0.0040*x39 - 0.2*u39 = 0;

solve;
display x0;
display y0;
display x1;
display y1;
display x2;
display y2;
display x3;
display y3;
display x4;
display y4;
display x5;
display y5;
display x6;
display y6;
display x7;
display y7;
display x8;
display y8;
display x9;
display y9;
display x10;
display y10;
display x11;
display y11;
display x12;
display y12;
display x13;
display y13;
display x14;
display y14;
display x15;
display y15;
display x16;
display y16;
display x17;
display y17;
display x18;
display y18;
display x19;
display y19;
display x20;
display y20;
display x21;
display y21;
display x22;
display y22;
display x23;
display y23;
display x24;
display y24;
display x25;
display y25;
display x26;
display y26;
display x27;
display y27;
display x28;
display y28;
display x29;
display y29;
display x30;
display y30;
display x31;
display y31;
display x32;
display y32;
display x33;
display y33;
display x34;
display y34;
display x35;
display y35;
display x36;
display y36;
display x37;
display y37;
display x38;
display y38;
display x39;
display y39;
display x40;
display y40;
display u0;
display u1;
display u2;
display u3;
display u4;
display u5;
display u6;
display u7;
display u8;
display u9;
display u10;
display u11;
display u12;
display u13;
display u14;
display u15;
display u16;
display u17;
display u18;
display u19;
display u20;
display u21;
display u22;
display u23;
display u24;
display u25;
display u26;
display u27;
display u28;
display u29;
display u30;
display u31;
display u32;
display u33;
display u34;
display u35;
display u36;
display u37;
display u38;
display u39;
display obj;
