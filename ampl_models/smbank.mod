#
#**************************
# SET UP THE INITIAL DATA *
#**************************
#   Problem :
#   *********
#   The small Bank Balancing problem (Thai model).
#   The problem is also named "MB64" in some references.
#   This is a nonlinear network problem with  conditioning
#   of the order of 10**4.
#   Source:
#   R. Dembo,
#   private communication, 1986.
#   SIF input: Ph. Toint, June 1990.
#   classification ONR2-MN-117-64
#   Number of arcs
#   Number of nodes
#   Constants
#   Objective
#   Network constraints
#   Network arcs
#   Solution
	param narcs := 117;
	param nodes := 64;

	var x1 >= 0.1 ,    := 0.1;
	var x2 >= 0.1 ,    := 0.1;
	var x3 >= 0.1 ,    := 0.1;
	var x4 >= 0.1 ,    := 0.7;
	var x5 >= 0.1 ,    := 0.7;
	var x6 >= 0.1 ,    := 0.1;
	var x7 >= 0.1 ,    := 0.1;
	var x8 >= 0.1 ,    := 0.1;
	var x9 >= 0.1 ,    := 0.1;
	var x10 >= 0.1 ,    := 0.1;
	var x11 >= 0.1 ,    := 0.1;
	var x12 >= 0.1 ,    := 0.7;
	var x13 >= 0.1 ,    := 0.1;
	var x14 >= 0.1 ,    := 0.1;
	var x15 >= 0.1 ,    := 0.1;
	var x16 >= 0.1 ,    := 0.1;
	var x17 >= 0.1 ,    := 0.9;
	var x18 >= 0.1 ,    := 0.1;
	var x19 >= 0.1 ,    := 0.1;
	var x20 >= 0.1 ,    := 0.1;
	var x21 >= 0.1 ,    := 0.6;
	var x22 >= 0.1 ,    := 0.9;
	var x23 >= 0.1 ,    := 0.1;
	var x24 >= 0.1 ,    := 0.1;
	var x25 >= 0.1 ,    := 0.1;
	var x26 >= 0.1 ,    := 0.1;
	var x27 >= 0.1 ,    := 0.1;
	var x28 >= 0.1 ,    := 0.4;
	var x29 >= 0.1 ,    := 0.1;
	var x30 >= 0.1 ,    := 0.4;
	var x31 >= 0.1 ,    := 0.4;
	var x32 >= 0.1 ,    := 0.1;
	var x33 >= 0.1 ,    := 1.0;
	var x34 >= 0.1 ,    := 0.5;
	var x35 >= 0.1 ,    := 0.5;
	var x36 >= 0.1 ,    := 0.5;
	var x37 >= 0.1 ,    := 1.1;
	var x38 >= 0.1 ,    := 0.1;
	var x39 >= 0.1 ,    := 0.1;
	var x40 >= 0.1 ,    := 0.1;
	var x41 >= 0.1 ,    := 0.1;
	var x42 >= 0.1 ,    := 0.3;
	var x43 >= 0.1 ,    := 0.1;
	var x44 >= 0.1 ,    := 0.1;
	var x45 >= 0.1 ,    := 0.1;
	var x46 >= 0.1 ,    := 0.1;
	var x47 >= 0.1 ,    := 0.2;
	var x48 >= 0.1 ,    := 0.1;
	var x49 >= 0.1 ,    := 0.1;
	var x50 >= 0.1 ,    := 0.1;
	var x51 >= 0.1 ,    := 0.1;
	var x52 >= 0.1 ,    := 0.1;
	var x53 >= 0.1 ,    := 0.1;
	var x54 >= 0.1 ,    := 0.1;
	var x55 >= 0.1 ,    := 0.1;
	var x56 >= 0.1 ,    := 0.1;
	var x57 >= 0.1 ,    := 0.1;
	var x58 >= 0.1 ,    := 0.1;
	var x59 >= 0.1 ,    := 0.1;
	var x60 >= 0.1 ,    := 0.1;
	var x61 >= 0.1 ,    := 0.7;
	var x62 >= 0.1 ,    := 0.1;
	var x63 >= 0.1 ,    := 0.2;
	var x64 >= 0.1 ,    := 0.3;
	var x65 >= 0.1 ,    := 0.2;
	var x66 >= 0.1 ,    := 0.1;
	var x67 >= 0.1 ,    := 0.2;
	var x68 >= 0.1 ,    := 0.2;
	var x69 >= 0.1 ,    := 0.1;
	var x70 >= 0.1 ,    := 0.1;
	var x71 >= 0.1 ,    := 0.1;
	var x72 >= 0.1 ,    := 0.1;
	var x73 >= 0.1 ,    := 0.1;
	var x74 >= 0.1 ,    := 0.1;
	var x75 >= 0.1 ,    := 0.1;
	var x76 >= 0.1 ,    := 0.1;
	var x77 >= 0.1 ,    := 0.1;
	var x78 >= 0.1 ,    := 0.1;
	var x79 >= 0.1 ,    := 0.1;
	var x80 >= 0.1 ,    := 0.1;
	var x81 >= 0.1 ,    := 0.2;
	var x82 >= 0.1 ,    := 0.1;
	var x83 >= 0.1 ,    := 0.1;
	var x84 >= 0.1 ,    := 0.1;
	var x85 ,    := 0.1;
	var x86 ,    := 0.1;
	var x87 ,    := 0.1;
	var x88 ,    := 0.7;
	var x89 ,    := 0.9;
	var x90 ,    := 1.0;
	var x91 ,    := 0.3;
	var x92 ,    := 1.3;
	var x93 ,    := 0.6;
	var x94 ,    := 0.9;
	var x95 ,    := 0.4;
	var x96 ,    := 0.5;
	var x97 ,    := 0.5;
	var x98 ,    := 0.4;
	var x99 ,    := 1.1;
	var x100 ,    := 0.5;
	var x101 ,    := 0.5;
	var x102 ,    := 0.5;
	var x103 ,    := 1.1;
	var x104 ,    := 0.7;
	var x105 ,    := 0.8;
	var x106 ,    := 0.7;
	var x107 ,    := 1.2;
	var x108 ,    := 0.2;
	var x109 ,    := 0.3;
	var x110 ,    := 0.2;
	var x111 ,    := 0.1;
	var x112 ,    := 0.2;
	var x113 ,    := 0.2;
	var x114 ,    := 0.1;
	var x115 ,    := 0.9;
	var x116 ,    := 0.7;
	var x117 >= 0.1 ,    := 0.1;

minimize obj:
	x1 * ((log(x1/175726.0)) - 1.0 )  + x2 * ((log(x2/140841.0)) - 1.0 )  + x3 * 
	((log(x3/12891.0)) - 1.0 )  + x4 * ((log(x4/273191.0)) - 1.0 )  + x5 * 
	((log(x5/273191.0)) - 1.0 )  + x6 * ((log(x6/12891.0)) - 1.0 )  + x7 * 
	((log(x7/140841.0)) - 1.0 )  + x8 * ((log(x8/175726.0)) - 1.0 )  + x9 * 
	((log(x9/4759.0)) - 1.0 )  + x10 * ((log(x10/380.0)) - 1.0 )  + x11 * 
	((log(x11/3859.0)) - 1.0 )  + x12 * ((log(x12/520758.0)) - 1.0 )  + x13 * 
	((log(x13/9940.0)) - 1.0 )  + x14 * ((log(x14/6164.0)) - 1.0 )  + x15 * 
	((log(x15/63118.0)) - 1.0 )  + x16 * ((log(x16/2738.0)) - 1.0 )  + x17 * 
	((log(x17/71083.0)) - 1.0 )  + x18 * ((log(x18/10153)) - 1.0 )  + 
	x19 * ((log(x19/9294.0)) - 1.0 )  + x20 * ((log(x20/4367.0)) - 1.0 )  + x21 * 
	((log(x21/442511.0)) - 1.0 )  + x22 * ((log(x22/84668)) - 1.0 )  + 
	x23 * ((log(x23/48524.0)) - 1.0 )  + x24 * ((log(x24/2645.0)) - 1.0 )  + x25 * 
	((log(x25/65210.0)) - 1.0 )  + x26 * ((log(x26/71528.0)) - 1.0 )  + x27 * 
	((log(x27/74346.0)) - 1.0 )  + x28 * ((log(x28/226846.0)) - 1.0 )  + x29 * 
	((log(x29/58725)) - 1.0 )  + x30 * ((log(x30/410017)) - 
	1.0 )  + x31 * ((log(x31/52774.0)) - 1.0 )  + x32 * ((log(x32/31847.0)) - 1.0 ) 
	 + x33 * ((log(x33/416672.0)) - 1.0 )  + x34 * ((log(x34/228906.0)) - 1.0 )  + 
	x35 * ((log(x35/437902)) - 1.0 )  + x36 * ((log(x36/53520.0)) - 1.0 
	)  + x37 * ((log(x37/434014.0)) - 1.0 )  + x38 * ((log(x38/11928.0)) - 1.0 )  + 
	x39 * ((log(x39/46747)) - 1.0 )  + x40 * ((log(x40/22284.0)) - 1.0 
	)  + x41 * ((log(x41/19093.0)) - 1.0 )  + x42 * ((log(x42/131949.0)) - 1.0 )  + 
	x43 * ((log(x43/29947.0)) - 1.0 )  + x44 * ((log(x44/982.0)) - 1.0 )  + x45 * 
	((log(x45/175800.0)) - 1.0 )  + x46 * ((log(x46/30274.0)) - 1.0 )  + x47 * 
	((log(x47/161094.0)) - 1.0 )  + x48 * ((log(x48/7203.0)) - 1.0 )  + x49 * 
	((log(x49/173359.0)) - 1.0 )  + x50 * ((log(x50/18971.0)) - 1.0 )  + x51 * 
	((log(x51/32905.0)) - 1.0 )  + x52 * ((log(x52/21978.0)) - 1.0 )  + x53 * 
	((log(x53/9613.0)) - 1.0 )  + x54 * ((log(x54/7720.0)) - 1.0 )  + x55 * 
	((log(x55/2846.0)) - 1.0 )  + x56 * ((log(x56/19522.0)) - 1.0 )  + x57 * 
	((log(x57/114482.0)) - 1.0 )  + x58 * ((log(x58/5996.0)) - 1.0 )  + x59 * 
	((log(x59/83376.0)) - 1.0 )  + x60 * ((log(x60/63295.0)) - 1.0 )  + x61 * 
	((log(x61/74619.0)) - 1.0 )  + x62 * ((log(x62/117681)) - 1.0 )  + 
	x63 * ((log(x63/3095.0)) - 1.0 )  + x64 * ((log(x64/140757.0)) - 1.0 )  + x65 * 
	((log(x65/60035.0)) - 1.0 )  + x66 * ((log(x66/25435.0)) - 1.0 )  + x67 * 
	((log(x67/77587.0)) - 1.0 )  + x68 * ((log(x68/58884.0)) - 1.0 )  + x69 * 
	((log(x69/31847.0)) - 1.0 )  + x70 * ((log(x70/159.0)) - 1.0 )  + x71 * 
	((log(x71/3241.0)) - 1.0 )  + x72 * ((log(x72/2138.0)) - 1.0 )  + x73 * 
	((log(x73/16583.0)) - 1.0 )  + x74 * ((log(x74/929)) - 1.0 )  + 
	x75 * ((log(x75/17342.0)) - 1.0 )  + x76 * ((log(x76/746.0)) - 1.0 )  + x77 * 
	((log(x77/27885.0)) - 1.0 )  + x78 * ((log(x78/2060.0)) - 1.0 )  + x79 * 
	((log(x79/25435.0)) - 1.0 )  + x80 * ((log(x80/57897.0)) - 1.0 )  + x81 * 
	((log(x81/124174.0)) - 1.0 )  + x82 * ((log(x82/2166.0)) - 1.0 )  + x83 * 
	((log(x83/2.0)) - 1.0 )  + x84 * ((log(x84/259.0)) - 1.0 )  + x117 * 
	((log(x117/14406)) - 1.0 ) ;

subject to n1:
	-x1 + x85 = 0;
subject to n2:
	-x2 + x86 = 0;
subject to n3:
	-x3 + x87 = 0;
subject to n4:
	-x4 + x88 = 0;
subject to n5:
	-x5 - x6 - x7 - x8 + x89 + x117 = 0;
subject to n6:
	-x9 - x10 - x11 - x12 + x90 = 0;
subject to n7:
	-x13 - x14 - x15 + x91 = 0;
subject to n8:
	-x16 - x17 - x18 - x19 - x20 + x92 = 0;
subject to n9:
	-x21 + x93 = 0;
subject to n10:
	-x22 + x94 = 0;
subject to n11:
	-x23 - x24 - x25 - x26 + x95 = 0;
subject to n12:
	-x27 - x28 + x96 = 0;
subject to n13:
	-x29 - x30 + x97 = 0;
subject to n14:
	-x31 + x98 = 0;
subject to n15:
	-x32 - x33 + x99 = 0;
subject to n16:
	-x34 + x100 = 0;
subject to n17:
	-x35 + x101 = 0;
subject to n18:
	-x36 + x102 = 0;
subject to n19:
	-x37 + x103 = 0;
subject to n20:
	-x38 - x39 - x40 - x41 - x42 + x104 = 0;
subject to n21:
	-x43 - x44 - x45 - x46 - x47 - x48 - x49 + x105 = 0;
subject to n22:
	-x50 - x51 - x52 - x53 - x54 - x55 - x56 + x106 = 0;
subject to n23:
	-x57 - x58 - x59 - x60 - x61 - x62 + x107 = 0;
subject to n24:
	-x63 + x108 = 0;
subject to n25:
	-x64 + x109 = 0;
subject to n26:
	-x65 + x110 = 0;
subject to n27:
	-x66 + x111 = 0;
subject to n28:
	-x67 + x112 = 0;
subject to n29:
	-x68 + x113 = 0;
subject to n30:
	-x69 + x114 = 0;
subject to n31:
	-x70 - x71 - x72 - x73 - x74 - x75 - x76 - x77 - x78 + x115 = 0;
subject to n32:
	-x79 - x80 - x81 - x82 - x83 - x84 + x116 = 0;
subject to n33:
	x8 - x85 = 0;
subject to n34:
	x7 - x86 = 0;
subject to n35:
	x6 - x87 = 0;
subject to n36:
	x5 - x88 = 0;
subject to n37:
	x12 + x15 + x20 - x89 = 0;
subject to n38:
	x14 + x19 + x21 + x26 + x84 - x90 = 0;
subject to n39:
	x11 + x18 + x25 - x91 = 0;
subject to n40:
	x10 + x13 + x22 + x24 + x83 - x92 = 0;
subject to n41:
	x42 + x49 + x56 + x62 - x93 = 0;
subject to n42:
	x48 + x55 + x61 - x94 = 0;
subject to n43:
	x41 + x47 + x54 - x95 = 0;
subject to n44:
	x1 + x40 + x46 + x53 + x60 - x96 = 0;
subject to n45:
	x2 + x39 + x45 + x52 + x59 - x97 = 0;
subject to n46:
	x3 + x44 + x51 + x58 - x98 = 0;
subject to n47:
	x4 + x38 + x43 + x50 + x57 - x99 = 0;
subject to n48:
	x28 + x78 - x100 = 0;
subject to n49:
	x30 + x77 - x101 = 0;
subject to n50:
	x31 + x76 - x102 = 0;
subject to n51:
	x33 + x75 - x103 = 0;
subject to n52:
	x34 + x63 - x104 = 0;
subject to n53:
	x35 + x64 - x105 = 0;
subject to n54:
	x36 + x65 - x106 = 0;
subject to n55:
	x37 + x66 - x107 = 0;
subject to n56:
	x74 + x82 - x108 = 0;
subject to n57:
	x73 + x81 - x109 = 0;
subject to n58:
	x72 + x80 - x110 = 0;
subject to n59:
	x79 - x111 = 0;
subject to n60:
	x27 + x71 - x112 = 0;
subject to n61:
	x29 + x70 - x113 = 0;
subject to n62:
	x32 - x114 = 0;
subject to n63:
	x17 - x115 = 0;
subject to n64:
	x9 + x16 + x23 + x67 + x68 + x69 - x116 - x117 = 0;

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
	display x61;
	display x62;
	display x63;
	display x64;
	display x65;
	display x66;
	display x67;
	display x68;
	display x69;
	display x70;
	display x71;
	display x72;
	display x73;
	display x74;
	display x75;
	display x76;
	display x77;
	display x78;
	display x79;
	display x80;
	display x81;
	display x82;
	display x83;
	display x84;
	display x85;
	display x86;
	display x87;
	display x88;
	display x89;
	display x90;
	display x91;
	display x92;
	display x93;
	display x94;
	display x95;
	display x96;
	display x97;
	display x98;
	display x99;
	display x100;
	display x101;
	display x102;
	display x103;
	display x104;
	display x105;
	display x106;
	display x107;
	display x108;
	display x109;
	display x110;
	display x111;
	display x112;
	display x113;
	display x114;
	display x115;
	display x116;
	display x117;
display obj;
