#
#**************************
# SET UP THE INITIAL DATA *
#**************************
#   Problem :
#   *********
#   A nonlinear minmax problem.
#   Source:
#   K. Jonasson and K. Madsen,
#   "Corrected sequential linear programming for sparse
#   minimax optimization", Technical report, Institute for Numerical
#   Analysis, Technical U. of Denmark.
#   SIF input: Nick Gould, October 1992.
#   classification LOR2-AN-V-V
#  the number of functions
#IE M                   3
#IE M                   8
#IE M                   14
#  the number of independent variables
#  other useful parameters
#  Define constants
#   The independent variables
param m := 20;
param n := 3 * (20);
param nm3 := -3 + (3 * (20));
param nm5 := -5 + (3 * (20));
param im5 := -5 + (57);
param ip3 := 3 + (57);
param id3 := 3;
param im1 := -1 + (60);
param im2 := -2 + (60);

var x1;
var x2;
var x3;
var x4;
var x5;
var x6;
var x7;
var x8;
var x9;
var x10;
var x11;
var x12;
var x13;
var x14;
var x15;
var x16;
var x17;
var x18;
var x19;
var x20;
var x21;
var x22;
var x23;
var x24;
var x25;
var x26;
var x27;
var x28;
var x29;
var x30;
var x31;
var x32;
var x33;
var x34;
var x35;
var x36;
var x37;
var x38;
var x39;
var x40;
var x41;
var x42;
var x43;
var x44;
var x45;
var x46;
var x47;
var x48;
var x49;
var x50;
var x51;
var x52;
var x53;
var x54;
var x55;
var x56;
var x57;
var x58;
var x59;
var x60;
var f;

minimize obj:
	f;

subject to c1:
	0 >= x3 * x3 + (cosh(x2)) + 2.0d+0 * x1 * x1 * x3 - f - 2.0*x3 - x6;
subject to c2:
	0 >= x6 * x6 + (cosh(x5)) + 2.0d+0 * x4 * x4 * x6 - f + x1 - 2.0*x6 - x9;
subject to c3:
	0 >= x9 * x9 + (cosh(x8)) + 2.0d+0 * x7 * x7 * x9 - f + x4 - 2.0*x9 - x12;
subject to c4:
	0 >= x12 * x12 + (cosh(x11)) + 2.0d+0 * x10 * x10 * x12 - f + x7 - 2.0*x12 - x15;
subject to c5:
	0 >= x15 * x15 + (cosh(x14)) + 2.0d+0 * x13 * x13 * x15 - f + x10 - 2.0*x15 - x18;
subject to c6:
	0 >= x18 * x18 + (cosh(x17)) + 2.0d+0 * x16 * x16 * x18 - f + x13 - 2.0*x18 - x21;
subject to c7:
	0 >= x21 * x21 + (cosh(x20)) + 2.0d+0 * x19 * x19 * x21 - f + x16 - 2.0*x21 - x24;
subject to c8:
	0 >= x24 * x24 + (cosh(x23)) + 2.0d+0 * x22 * x22 * x24 - f + x19 - 2.0*x24 - x27;
subject to c9:
	0 >= x27 * x27 + (cosh(x26)) + 2.0d+0 * x25 * x25 * x27 - f + x22 - 2.0*x27 - x30;
subject to c10:
	0 >= x30 * x30 + (cosh(x29)) + 2.0d+0 * x28 * x28 * x30 - f + x25 - 2.0*x30 - x33;
subject to c11:
	0 >= x33 * x33 + (cosh(x32)) + 2.0d+0 * x31 * x31 * x33 - f + x28 - 2.0*x33 - x36;
subject to c12:
	0 >= x36 * x36 + (cosh(x35)) + 2.0d+0 * x34 * x34 * x36 - f + x31 - 2.0*x36 - x39;
subject to c13:
	0 >= x39 * x39 + (cosh(x38)) + 2.0d+0 * x37 * x37 * x39 - f + x34 - 2.0*x39 - x42;
subject to c14:
	0 >= x42 * x42 + (cosh(x41)) + 2.0d+0 * x40 * x40 * x42 - f + x37 - 2.0*x42 - x45;
subject to c15:
	0 >= x45 * x45 + (cosh(x44)) + 2.0d+0 * x43 * x43 * x45 - f + x40 - 2.0*x45 - x48;
subject to c16:
	0 >= x48 * x48 + (cosh(x47)) + 2.0d+0 * x46 * x46 * x48 - f + x43 - 2.0*x48 - x51;
subject to c17:
	0 >= x51 * x51 + (cosh(x50)) + 2.0d+0 * x49 * x49 * x51 - f + x46 - 2.0*x51 - x54;
subject to c18:
	0 >= x54 * x54 + (cosh(x53)) + 2.0d+0 * x52 * x52 * x54 - f + x49 - 2.0*x54 - x57;
subject to c19:
	0 >= x57 * x57 + (cosh(x56)) + 2.0d+0 * x55 * x55 * x57 - f + x52 - 2.0*x57 - x60;
subject to c20:
	0 >= x60 * x60 + (cosh(x59)) + 2.0d+0 * x58 * x58 * x60 - f + x55 - 2.0*x60;

solve;
display x1;
display x2;
display x3;
display x4;
display x5;
display x6;
display x7;
display x8;
display x9;
display x10;
display x11;
display x12;
display x13;
display x14;
display x15;
display x16;
display x17;
display x18;
display x19;
display x20;
display x21;
display x22;
display x23;
display x24;
display x25;
display x26;
display x27;
display x28;
display x29;
display x30;
display x31;
display x32;
display x33;
display x34;
display x35;
display x36;
display x37;
display x38;
display x39;
display x40;
display x41;
display x42;
display x43;
display x44;
display x45;
display x46;
display x47;
display x48;
display x49;
display x50;
display x51;
display x52;
display x53;
display x54;
display x55;
display x56;
display x57;
display x58;
display x59;
display x60;
display f;
display obj;
