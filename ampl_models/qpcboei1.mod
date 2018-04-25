# AMPL Model by Hande Y. Benson
#
# Copyright (C) 2001 Princeton University
# All Rights Reserved
#
# Permission to use, copy, modify, and distribute this software and
# its documentation for any purpose and without fee is hereby
# granted, provided that the above copyright notice appear in all
# copies and that the copyright notice and this
# permission notice appear in all supporting documentation.                     

#   Source: a variant on the BOEING1 linear programming problem
#   with an additional CONVEX diagonal Hessian matrix as given by
#   N. I. M. Gould, "An algorithm for large-scale quadratic programming",
#   IMA J. Num. Anal (1991), 11, 299-324, problem class 4.

#   SIF input: Nick Gould, January 1993

#   classification QLR2-MN-384-351

#
#   Problem :
#   *********
#   Source: a variant on the BOEING1 linear programming problem
#   with an additional CONVEX diagonal Hessian matrix as given by
#   N. I. M. Gould, "An algorithm for large-scale quadratic programming",
#   IMA J. Num. Anal (1991), 11, 299-324, problem class 4.
#   SIF input: Nick Gould, January 1993
	param n := 384;

	var pboshnl0 >= 0.0, := 0.0;
	var pboshnl1 >= 0.0, := 0.0;
	var pboshnl7 >= 0.0, := 0.0;
	var pboshnl8 >= 0.0, := 0.0;
	var pboslax0 >= 0.0, := 0.0;
	var pboslax1 >= 0.0, := 0.0;
	var pboslax7 >= 0.0, := 0.0;
	var pbossea0 >= 0.0, := 0.0;
	var pbossea1 >= 0.0, := 0.0;
	var pbossea2 >= 0.0, := 0.0;
	var pbossfo0 >= 0.0, := 0.0;
	var pbossfo1 >= 0.0, := 0.0;
	var pbostpe1 >= 0.0, := 0.0;
	var pbostpe2 >= 0.0, := 0.0;
	var pbostyo1 >= 0.0, := 0.0;
	var pbostyo2 >= 0.0, := 0.0;
	var pbosyul0 >= 0.0, := 0.0;
	var pbosyul1 >= 0.0, := 0.0;
	var pbosyul2 >= 0.0, := 0.0;
	var pbosyul3 >= 0.0, := 0.0;
	var pbosyul4 >= 0.0, := 0.0;
	var pbosyvr0 >= 0.0, := 0.0;
	var pbosyvr1 >= 0.0, := 0.0;
	var pbosyvr2 >= 0.0, := 0.0;
	var pbosywg0 >= 0.0, := 0.0;
	var pbosywg1 >= 0.0, := 0.0;
	var pbosywg7 >= 0.0, := 0.0;
	var pbosyyz0 >= 0.0, := 0.0;
	var pbosyyz1 >= 0.0, := 0.0;
	var pburoak0 >= 0.0, := 0.0;
	var pburoak1 >= 0.0, := 0.0;
	var pbursea0 >= 0.0, := 0.0;
	var pbursea1 >= 0.0, := 0.0;
	var pbursfo0 >= 0.0, := 0.0;
	var phnllax0 >= 0.0, := 0.0;
	var phnllax1 >= 0.0, := 0.0;
	var phnllax2 >= 0.0, := 0.0;
	var phnllax3 >= 0.0, := 0.0;
	var phnllon0 >= 0.0, := 0.0;
	var phnllon6 >= 0.0, := 0.0;
	var phnlpar0 >= 0.0, := 0.0;
	var phnlpar6 >= 0.0, := 0.0;
	var phnlsea0 >= 0.0, := 0.0;
	var phnlsea1 >= 0.0, := 0.0;
	var phnlsea2 >= 0.0, := 0.0;
	var phnlsfo0 >= 0.0, := 0.0;
	var phnlsfo1 >= 0.0, := 0.0;
	var phnlsfo7 >= 0.0, := 0.0;
	var phnlyvr0 >= 0.0, := 0.0;
	var phnlyvr1 >= 0.0, := 0.0;
	var phnlyvr7 >= 0.0, := 0.0;
	var plassea0 >= 0.0, := 0.0;
	var plassea1 >= 0.0, := 0.0;
	var plasyvr0 >= 0.0, := 0.0;
	var plasyvr6 >= 0.0, := 0.0;
	var plaxoak0 >= 0.0, := 0.0;
	var plaxoak1 >= 0.0, := 0.0;
	var plaxoak2 >= 0.0, := 0.0;
	var plaxsea0 >= 0.0, := 0.0;
	var plaxsea1 >= 0.0, := 0.0;
	var plaxsea2 >= 0.0, := 0.0;
	var plaxsea3 >= 0.0, := 0.0;
	var plaxsea4 >= 0.0, := 0.0;
	var plaxsea5 >= 0.0, := 0.0;
	var plaxsea6 >= 0.0, := 0.0;
	var plaxsea7 >= 0.0, := 0.0;
	var plaxsea8 >= 0.0, := 0.0;
	var plaxsea9 >= 0.0, := 0.0;
	var plaxsfo0 >= 0.0, := 0.0;
	var plaxsfo1 >= 0.0, := 0.0;
	var plaxsfo2 >= 0.0, := 0.0;
	var plaxsfo3 >= 0.0, := 0.0;
	var plaxsfo4 >= 0.0, := 0.0;
	var plaxsfo5 >= 0.0, := 0.0;
	var plaxsfo6 >= 0.0, := 0.0;
	var plaxtpe0 >= 0.0, := 0.0;
	var plaxtpe1 >= 0.0, := 0.0;
	var plaxtpe2 >= 0.0, := 0.0;
	var plaxtpe8 >= 0.0, := 0.0;
	var plaxtyo0 >= 0.0, := 0.0;
	var plaxtyo1 >= 0.0, := 0.0;
	var plaxtyo2 >= 0.0, := 0.0;
	var plaxtyo8 >= 0.0, := 0.0;
	var plaxyvr0 >= 0.0, := 0.0;
	var plaxyvr1 >= 0.0, := 0.0;
	var plaxyvr2 >= 0.0, := 0.0;
	var plonpar0 >= 0.0, := 0.0;
	var plonpar1 >= 0.0, := 0.0;
	var plonpar2 >= 0.0, := 0.0;
	var plonsea0 >= 0.0, := 0.0;
	var plonsea1 >= 0.0, := 0.0;
	var plonyvr0 >= 0.0, := 0.0;
	var plonyvr1 >= 0.0, := 0.0;
	var plonyvr7 >= 0.0, := 0.0;
	var poakont0 >= 0.0, := 0.0;
	var poakont1 >= 0.0, := 0.0;
	var poaksea0 >= 0.0, := 0.0;
	var poaksea1 >= 0.0, := 0.0;
	var poaksea2 >= 0.0, := 0.0;
	var pontsfo0 >= 0.0, := 0.0;
	var pontsea0 >= 0.0, := 0.0;
	var pontsea1 >= 0.0, := 0.0;
	var pontsea2 >= 0.0, := 0.0;
	var pparsea0 >= 0.0, := 0.0;
	var pparsea1 >= 0.0, := 0.0;
	var pparyvr0 >= 0.0, := 0.0;
	var pparyvr1 >= 0.0, := 0.0;
	var pparyvr7 >= 0.0, := 0.0;
	var prnosea0 >= 0.0, := 0.0;
	var prnosea1 >= 0.0, := 0.0;
	var prnoyvr0 >= 0.0, := 0.0;
	var pseasfo0 >= 0.0, := 0.0;
	var pseasfo1 >= 0.0, := 0.0;
	var pseasfo2 >= 0.0, := 0.0;
	var pseasfo3 >= 0.0, := 0.0;
	var pseasfo4 >= 0.0, := 0.0;
	var pseasfo5 >= 0.0, := 0.0;
	var pseasfo6 >= 0.0, := 0.0;
	var pseatpe0 >= 0.0, := 0.0;
	var pseatpe1 >= 0.0, := 0.0;
	var pseatpe2 >= 0.0, := 0.0;
	var pseatpe3 >= 0.0, := 0.0;
	var pseatyo0 >= 0.0, := 0.0;
	var pseatyo1 >= 0.0, := 0.0;
	var pseatyo2 >= 0.0, := 0.0;
	var pseatyo3 >= 0.0, := 0.0;
	var pseayvr0 >= 0.0, := 0.0;
	var pseayvr1 >= 0.0, := 0.0;
	var pseayvr2 >= 0.0, := 0.0;
	var pseayvr3 >= 0.0, := 0.0;
	var pseayvr4 >= 0.0, := 0.0;
	var pseayvr5 >= 0.0, := 0.0;
	var pseayvr6 >= 0.0, := 0.0;
	var pseayvr7 >= 0.0, := 0.0;
	var pseayvr8 >= 0.0, := 0.0;
	var psfotpe0 >= 0.0, := 0.0;
	var psfotpe1 >= 0.0, := 0.0;
	var psfotpe2 >= 0.0, := 0.0;
	var psfotpe8 >= 0.0, := 0.0;
	var psfotyo0 >= 0.0, := 0.0;
	var psfotyo1 >= 0.0, := 0.0;
	var psfotyo2 >= 0.0, := 0.0;
	var psfotyo8 >= 0.0, := 0.0;
	var psfoyvr0 >= 0.0, := 0.0;
	var psfoyvr1 >= 0.0, := 0.0;
	var ptpetyo0 >= 0.0, := 0.0;
	var ptpetyo1 >= 0.0, := 0.0;
	var ptpetyo2 >= 0.0, := 0.0;
	var ptpetyo3 >= 0.0, := 0.0;
	var ptpeyvr0 >= 0.0, := 0.0;
	var ptyoyvr0 >= 0.0, := 0.0;
	var pyulyvr0 >= 0.0, := 0.0;
	var pyulyvr1 >= 0.0, := 0.0;
	var pyulyvr2 >= 0.0, := 0.0;
	var pyulyvr3 >= 0.0, := 0.0;
	var pyulywg0 >= 0.0, := 0.0;
	var pyulywg1 >= 0.0, := 0.0;
	var pyulywg2 >= 0.0, := 0.0;
	var pyulywg3 >= 0.0, := 0.0;
	var pyulyyz0 >= 0.0, := 0.0;
	var pyulyyz1 >= 0.0, := 0.0;
	var pyulyyz2 >= 0.0, := 0.0;
	var pyulyyz3 >= 0.0, := 0.0;
	var pyulyyz4 >= 0.0, := 0.0;
	var pyvrywg0 >= 0.0, := 0.0;
	var pyvrywg1 >= 0.0, := 0.0;
	var pyvrywg2 >= 0.0, := 0.0;
	var pyvryyz0 >= 0.0, := 0.0;
	var pyvryyz1 >= 0.0, := 0.0;
	var pyvryyz2 >= 0.0, := 0.0;
	var pywgyyz0 >= 0.0, := 0.0;
	var pywgyyz1 >= 0.0, := 0.0;
	var pywgyyz2 >= 0.0, := 0.0;
	var pywgyyz3 >= 0.0, := 0.0;
	var pbosoak0 >= 0.0, := 0.0;
	var pbosoak6 >= 0.0, := 0.0;
	var pbosbur1 >= 0.0, := 0.0;
	var pbosbur2 >= 0.0, := 0.0;
	var pbosont1 >= 0.0, := 0.0;
	var pbosont2 >= 0.0, := 0.0;
	var pburyvr1 >= 0.0, := 0.0;
	var pburtyo1 >= 0.0, := 0.0;
	var pburtpe1 >= 0.0, := 0.0;
	var pburhnl0 >= 0.0, := 0.0;
	var pburhnl6 >= 0.0, := 0.0;
	var phnloak0 >= 0.0, := 0.0;
	var phnloak1 >= 0.0, := 0.0;
	var phnloak2 >= 0.0, := 0.0;
	var phnloak8 >= 0.0, := 0.0;
	var phnlont0 >= 0.0, := 0.0;
	var phnlont6 >= 0.0, := 0.0;
	var phnlywg1 >= 0.0, := 0.0;
	var phnlyyz1 >= 0.0, := 0.0;
	var phnlyul1 >= 0.0, := 0.0;
	var plastyo1 >= 0.0, := 0.0;
	var plastpe1 >= 0.0, := 0.0;
	var plaxlon0 >= 0.0, := 0.0;
	var plaxlon6 >= 0.0, := 0.0;
	var plaxlon7 >= 0.0, := 0.0;
	var plaxpar0 >= 0.0, := 0.0;
	var plaxpar6 >= 0.0, := 0.0;
	var plaxpar7 >= 0.0, := 0.0;
	var pburlon1 >= 0.0, := 0.0;
	var pburpar1 >= 0.0, := 0.0;
	var plonont1 >= 0.0, := 0.0;
	var plonoak1 >= 0.0, := 0.0;
	var poakpar1 >= 0.0, := 0.0;
	var poaktyo1 >= 0.0, := 0.0;
	var poaktpe1 >= 0.0, := 0.0;
	var pontpar1 >= 0.0, := 0.0;
	var ponttyo1 >= 0.0, := 0.0;
	var ponttpe1 >= 0.0, := 0.0;
	var pparsfo1 >= 0.0, := 0.0;
	var prnotyo1 >= 0.0, := 0.0;
	var prnotpe1 >= 0.0, := 0.0;
	var ptpeywg1 >= 0.0, := 0.0;
	var ptpeyyz1 >= 0.0, := 0.0;
	var ptpeyul1 >= 0.0, := 0.0;
	var ptyoyul1 >= 0.0, := 0.0;
	var ptyoyyz1 >= 0.0, := 0.0;
	var ptyoywg1 >= 0.0, := 0.0;
	var plaxont0 >= 0.0, := 0.0;
	var grdtimo1 >= 0.0, := 0.0;
	var grdtimn1 >= -105.0 ,  <= 0.0;
	var grdtimo2 >= 0.0, := 0.0;
	var grdtimn2 >= -91.0 ,  <= 0.0;
	var grdtimo3 >= 0.0, := 0.0;
	var grdtimn3 >= -47.0 ,  <= 0.0;
	var grdtimo4 >= 0.0, := 0.0;
	var grdtimn4 >= -43.5 ,  <= 0.0;
	var grdtimo5 >= 0.0, := 0.0;
	var grdtimn5 >= -87.0 ,  <= 0.0;
	var grdtimo6 >= 0.0, := 0.0;
	var grdtimn6 >= -81.0 ,  <= 0.0;
	var n1001ac1 >= 0.0, := 0.0 ,  <= 3.0;
	var n1001ac2 >= 0.0, := 0.0 ,  <= 3.0;
	var n1001ac3 >= 0.0, := 0.0 ,  <= 3.0;
	var n1002ac1 >= 0.0, := 0.0 ,  <= 3.0;
	var n1002ac2 >= 0.0, := 0.0 ,  <= 3.0;
	var n1002ac3 >= 0.0, := 0.0 ,  <= 3.0;
	var n1003ac1 >= 0.0, := 0.0 ,  <= 4.0;
	var n1003ac2 >= 0.0, := 0.0 ,  <= 4.0;
	var n1003ac3 >= 0.0, := 0.0 ,  <= 4.0;
	var n1004ac1 >= 0.0, := 0.0 ,  <= 4.0;
	var n1004ac2 >= 0.0, := 0.0 ,  <= 4.0;
	var n1004ac3 >= 0.0, := 0.0 ,  <= 4.0;
	var n1005ac3 >= 0.0, := 0.0 ,  <= 2.0;
	var n1105ac3 >= 0.0, := 0.0 ,  <= 1.0;
	var n1006ac3 >= 0.0, := 0.0 ,  <= 2.0;
	var n1007ac1 >= 0.0, := 0.0 ,  <= 2.0;
	var n1007ac2 >= 0.0, := 0.0 ,  <= 2.0;
	var n1007ac3 >= 0.0, := 0.0 ,  <= 2.0;
	var n1008ac1 >= 0.0, := 0.0 ,  <= 7.0;
	var n1008ac2 >= 0.0, := 0.0 ,  <= 7.0;
	var n1008ac3 >= 0.0, := 0.0 ,  <= 7.0;
	var n1008ac4 >= 0.0, := 0.0 ,  <= 7.0;
	var n1008ac5 >= 0.0, := 0.0 ,  <= 7.0;
	var n1008ac6 >= 0.0, := 0.0 ,  <= 7.0;
	var n1009ac1 >= 0.0, := 0.0 ,  <= 7.0;
	var n1009ac2 >= 0.0, := 0.0 ,  <= 7.0;
	var n1009ac3 >= 0.0, := 0.0 ,  <= 7.0;
	var n1009ac4 >= 0.0, := 0.0 ,  <= 7.0;
	var n1009ac5 >= 0.0, := 0.0 ,  <= 7.0;
	var n1010ac1 >= 0.0, := 0.0 ,  <= 7.0;
	var n1010ac2 >= 0.0, := 0.0 ,  <= 7.0;
	var n1010ac3 >= 0.0, := 0.0 ,  <= 7.0;
	var n1010ac4 >= 0.0, := 0.0 ,  <= 7.0;
	var n1010ac5 >= 0.0, := 0.0 ,  <= 7.0;
	var n1010ac6 >= 0.0, := 0.0 ,  <= 7.0;
	var n1011ac1 >= 0.0, := 0.0 ,  <= 7.0;
	var n1011ac2 >= 0.0, := 0.0 ,  <= 7.0;
	var n1011ac3 >= 0.0, := 0.0 ,  <= 7.0;
	var n1011ac4 >= 0.0, := 0.0 ,  <= 7.0;
	var n1011ac5 >= 0.0, := 0.0 ,  <= 7.0;
	var n1011ac6 >= 0.0, := 0.0 ,  <= 7.0;
	var n1012ac1 >= 0.0, := 0.0 ,  <= 7.0;
	var n1012ac2 >= 0.0, := 0.0 ,  <= 7.0;
	var n1012ac3 >= 0.0, := 0.0 ,  <= 7.0;
	var n1012ac4 >= 0.0, := 0.0 ,  <= 7.0;
	var n1012ac5 >= 0.0, := 0.0 ,  <= 7.0;
	var n1013ac3 >= 0.0, := 0.0 ,  <= 4.0;
	var n1013ac4 >= 0.0, := 0.0 ,  <= 4.0;
	var n1013ac5 >= 0.0, := 0.0 ,  <= 4.0;
	var n1013ac6 >= 0.0, := 0.0 ,  <= 4.0;
	var n1014ac3 >= 0.0, := 0.0 ,  <= 4.0;
	var n1014ac4 >= 0.0, := 0.0 ,  <= 4.0;
	var n1014ac5 >= 0.0, := 0.0 ,  <= 4.0;
	var n1014ac6 >= 0.0, := 0.0 ,  <= 4.0;
	var n1015ac3 >= 0.0, := 0.0 ,  <= 4.0;
	var n1015ac4 >= 0.0, := 0.0 ,  <= 4.0;
	var n1015ac5 >= 0.0, := 0.0 ,  <= 4.0;
	var n1015ac6 >= 0.0, := 0.0 ,  <= 4.0;
	var n1016ac3 >= 0.0, := 0.0 ,  <= 4.0;
	var n1016ac4 >= 0.0, := 0.0 ,  <= 4.0;
	var n1016ac5 >= 0.0, := 0.0 ,  <= 4.0;
	var n1016ac6 >= 0.0, := 0.0 ,  <= 4.0;
	var n1017ac3 >= 0.0, := 0.0 ,  <= 4.0;
	var n1017ac4 >= 0.0, := 0.0 ,  <= 4.0;
	var n1017ac5 >= 0.0, := 0.0 ,  <= 4.0;
	var n1017ac6 >= 0.0, := 0.0 ,  <= 4.0;
	var n1018ac1 >= 0.0, := 0.0 ,  <= 7.0;
	var n1018ac2 >= 0.0, := 0.0 ,  <= 7.0;
	var n1018ac3 >= 0.0, := 0.0 ,  <= 7.0;
	var n1018ac4 >= 0.0, := 0.0 ,  <= 7.0;
	var n1018ac5 >= 0.0, := 0.0 ,  <= 7.0;
	var n1018ac6 >= 0.0, := 0.0 ,  <= 7.0;
	var n1019ac1 >= 0.0, := 0.0 ,  <= 7.0;
	var n1019ac2 >= 0.0, := 0.0 ,  <= 7.0;
	var n1019ac3 >= 0.0, := 0.0 ,  <= 7.0;
	var n1019ac4 >= 0.0, := 0.0 ,  <= 7.0;
	var n1019ac5 >= 0.0, := 0.0 ,  <= 7.0;
	var n1020ac1 >= 0.0, := 0.0 ,  <= 7.0;
	var n1020ac2 >= 0.0, := 0.0 ,  <= 7.0;
	var n1020ac3 >= 0.0, := 0.0 ,  <= 7.0;
	var n1020ac4 >= 0.0, := 0.0 ,  <= 7.0;
	var n1020ac5 >= 0.0, := 0.0 ,  <= 7.0;
	var n1020ac6 >= 0.0, := 0.0 ,  <= 7.0;
	var n1021ac1 >= 0.0, := 0.0 ,  <= 7.0;
	var n1021ac2 >= 0.0, := 0.0 ,  <= 7.0;
	var n1021ac3 >= 0.0, := 0.0 ,  <= 7.0;
	var n1021ac4 >= 0.0, := 0.0 ,  <= 7.0;
	var n1021ac5 >= 0.0, := 0.0 ,  <= 7.0;
	var n1022ac1 >= 0.0, := 0.0 ,  <= 1.0;
	var n1023ac1 >= 0.0, := 0.0 ,  <= 1.0;
	var n1026ac1 >= 0.0, := 0.0 ,  <= 1.0;
	var n1027ac1 >= 0.0, := 0.0 ,  <= 1.0;
	var n1028ac1 >= 0.0, := 0.0 ,  <= 1.0;
	var n1029ac1 >= 0.0, := 0.0 ,  <= 1.0;
	var n1030ac1 >= 0.0, := 0.0 ,  <= 1.0;
	var n1032ac1 >= 0.0, := 0.0 ,  <= 1.0;
	var n1032ac2 >= 0.0, := 0.0 ,  <= 1.0;
	var n1032ac3 >= 0.0, := 0.0 ,  <= 1.0;
	var n1032ac4 >= 0.0, := 0.0 ,  <= 1.0;
	var n1032ac5 >= 0.0, := 0.0 ,  <= 1.0;
	var n1033ac1 >= 0.0, := 0.0 ,  <= 5.0;
	var n1033ac2 >= 0.0, := 0.0 ,  <= 5.0;
	var n1033ac3 >= 0.0, := 0.0 ,  <= 5.0;
	var n1033ac4 >= 0.0, := 0.0 ,  <= 5.0;
	var n1033ac5 >= 0.0, := 0.0 ,  <= 5.0;
	var n1034ac1 >= 0.0, := 0.0 ,  <= 5.0;
	var n1034ac2 >= 0.0, := 0.0 ,  <= 5.0;
	var n1034ac3 >= 0.0, := 0.0 ,  <= 5.0;
	var n1035ac1 >= 0.0, := 0.0 ,  <= 5.0;
	var n1035ac2 >= 0.0, := 0.0 ,  <= 5.0;
	var n1035ac3 >= 0.0, := 0.0 ,  <= 5.0;
	var n1035ac4 >= 0.0, := 0.0 ,  <= 5.0;
	var n1035ac5 >= 0.0, := 0.0 ,  <= 5.0;
	var n1036ac1 >= 0.0, := 0.0 ,  <= 5.0;
	var n1036ac2 >= 0.0, := 0.0 ,  <= 5.0;
	var n1036ac3 >= 0.0, := 0.0 ,  <= 5.0;
	var n1037ac4 >= 0.0, := 0.0 ,  <= 5.0;
	var n1037ac5 >= 0.0, := 0.0 ,  <= 5.0;
	var n1038ac4 >= 0.0, := 0.0 ,  <= 10.0;
	var n1038ac5 >= 0.0, := 0.0 ,  <= 10.0;
	var n1039ac4 >= 0.0, := 0.0 ,  <= 7.0;
	var n1039ac5 >= 0.0, := 0.0 ,  <= 7.0;
	var n1040ac4 >= 0.0, := 0.0 ,  <= 10.0;
	var n1040ac5 >= 0.0, := 0.0 ,  <= 10.0;
	var n1040ac6 >= 0.0, := 0.0 ,  <= 10.0;
	var n1041ac4 >= 0.0, := 0.0 ,  <= 20.0;
	var n1041ac5 >= 0.0, := 0.0 ,  <= 20.0;
	var n1041ac6 >= 0.0, := 0.0 ,  <= 20.0;
	var n1042ac4 >= 0.0, := 0.0 ,  <= 20.0;
	var n1042ac5 >= 0.0, := 0.0 ,  <= 20.0;
	var n1042ac6 >= 0.0, := 0.0 ,  <= 20.0;
	var n1043ac1 >= 0.0, := 0.0 ,  <= 3.0;
	var n1043ac2 >= 0.0, := 0.0 ,  <= 3.0;
	var n1043ac3 >= 0.0, := 0.0 ,  <= 3.0;
	var n1044ac1 >= 0.0, := 0.0 ,  <= 3.0;
	var n1044ac2 >= 0.0, := 0.0 ,  <= 3.0;
	var n1044ac3 >= 0.0, := 0.0 ,  <= 3.0;
	var n1046ac3 >= 0.0, := 0.0 ,  <= 2.0;
	var n1047ac1 >= 0.0, := 0.0 ,  <= 1.0;
	var n1047ac2 >= 0.0, := 0.0 ,  <= 1.0;
	var n1047ac3 >= 0.0, := 0.0 ,  <= 1.0;
	var n1050ac3 >= 0.0, := 0.0 ,  <= 5.0;
	var n1050ac4 >= 0.0, := 0.0 ,  <= 5.0;
	var n1050ac5 >= 0.0, := 0.0 ,  <= 5.0;
	var n1051ac1 >= 0.0, := 0.0 ,  <= 20.0;
	var n1051ac2 >= 0.0, := 0.0 ,  <= 20.0;
	var n1051ac3 >= 0.0, := 0.0 ,  <= 20.0;
	var n1051ac4 >= 0.0, := 0.0 ,  <= 20.0;
	var n1051ac5 >= 0.0, := 0.0 ,  <= 20.0;
	var n1051ac6 >= 0.0, := 0.0 ,  <= 20.0;

minimize obj:
	1.0 * pboshnl0 * pboshnl0 + 1.023498654 * pboshnl1 * pboshnl1 + 1.046997428 * 
	pboshnl7 * pboshnl7 + 1.070496082 * pboshnl8 * pboshnl8 + 1.093994737 * 
	pboslax0 * pboslax0 + 1.11749351 * pboslax1 * pboslax1 + 1.140992165 * pboslax7 
	* pboslax7 + 1.164490819 * pbossea0 * pbossea0 + 1.187989593 * pbossea1 * 
	pbossea1 + 1.211488247 * pbossea2 * pbossea2 + 1.234986901 * pbossfo0 * 
	pbossfo0 + 1.258485675 * pbossfo1 * pbossfo1 + 1.281984329 * pbostpe1 * 
	pbostpe1 + 1.305482984 * pbostpe2 * pbostpe2 + 1.328981757 * pbostyo1 * 
	pbostyo1 + 1.352480412 * pbostyo2 * pbostyo2 + 1.375979066 * pbosyul0 * 
	pbosyul0 + 1.399477839 * pbosyul1 * pbosyul1 + 1.422976494 * pbosyul2 * 
	pbosyul2 + 1.446475148 * pbosyul3 * pbosyul3 + 1.469973922 * pbosyul4 * 
	pbosyul4 + 1.493472576 * pbosyvr0 * pbosyvr0 + 1.516971231 * pbosyvr1 * 
	pbosyvr1 + 1.540470004 * pbosyvr2 * pbosyvr2 + 1.563968658 * pbosywg0 * 
	pbosywg0 + 1.587467313 * pbosywg1 * pbosywg1 + 1.610966086 * pbosywg7 * 
	pbosywg7 + 1.634464741 * pbosyyz0 * pbosyyz0 + 1.657963395 * pbosyyz1 * 
	pbosyyz1 + 1.681462169 * pburoak0 * pburoak0 + 1.704960823 * pburoak1 * 
	pburoak1 + 1.728459477 * pbursea0 * pbursea0 + 1.751958251 * pbursea1 * 
	pbursea1 + 1.775456905 * pbursfo0 * pbursfo0 + 1.79895556 * phnllax0 * phnllax0 
	+ 1.822454333 * phnllax1 * phnllax1 + 1.845952988 * phnllax2 * phnllax2 + 
	1.869451642 * phnllax3 * phnllax3 + 1.892950416 * phnllon0 * phnllon0 + 
	1.91644907 * phnllon6 * phnllon6 + 1.939947724 * phnlpar0 * phnlpar0 + 
	1.963446498 * phnlpar6 * phnlpar6 + 1.986945152 * phnlsea0 * phnlsea0 + 
	2.010443926 * phnlsea1 * phnlsea1 + 2.033942461 * phnlsea2 * phnlsea2 + 
	2.057441235 * phnlsfo0 * phnlsfo0 + 2.080940008 * phnlsfo1 * phnlsfo1 + 
	2.104438543 * phnlsfo7 * phnlsfo7 + 2.127937317 * phnlyvr0 * phnlyvr0 + 
	2.15143609 * phnlyvr1 * phnlyvr1 + 2.174934626 * phnlyvr7 * phnlyvr7 + 
	2.198433399 * plassea0 * plassea0 + 2.221932173 * plassea1 * plassea1 + 
	2.245430708 * plasyvr0 * plasyvr0 + 2.268929482 * plasyvr6 * plasyvr6 + 
	2.292428255 * plaxoak0 * plaxoak0 + 2.31592679 * plaxoak1 * plaxoak1 + 
	2.339425564 * plaxoak2 * plaxoak2 + 2.362924337 * plaxsea0 * plaxsea0 + 
	2.386422873 * plaxsea1 * plaxsea1 + 2.409921646 * plaxsea2 * plaxsea2 + 
	2.43342042 * plaxsea3 * plaxsea3 + 2.456918955 * plaxsea4 * plaxsea4 + 
	2.480417728 * plaxsea5 * plaxsea5 + 2.503916502 * plaxsea6 * plaxsea6 + 
	2.527415037 * plaxsea7 * plaxsea7 + 2.550913811 * plaxsea8 * plaxsea8 + 
	2.574412584 * plaxsea9 * plaxsea9 + 2.597911119 * plaxsfo0 * plaxsfo0 + 
	2.621409893 * plaxsfo1 * plaxsfo1 + 2.644908667 * plaxsfo2 * plaxsfo2 + 
	2.668407202 * plaxsfo3 * plaxsfo3 + 2.691905975 * plaxsfo4 * plaxsfo4 + 
	2.715404749 * plaxsfo5 * plaxsfo5 + 2.738903284 * plaxsfo6 * plaxsfo6 + 
	2.762402058 * plaxtpe0 * plaxtpe0 + 2.785900831 * plaxtpe1 * plaxtpe1 + 
	2.809399366 * plaxtpe2 * plaxtpe2 + 2.83289814 * plaxtpe8 * plaxtpe8 + 
	2.856396914 * plaxtyo0 * plaxtyo0 + 2.879895449 * plaxtyo1 * plaxtyo1 + 
	2.903394222 * plaxtyo2 * plaxtyo2 + 2.926892996 * plaxtyo8 * plaxtyo8 + 
	2.950391531 * plaxyvr0 * plaxyvr0 + 2.973890305 * plaxyvr1 * plaxyvr1 + 
	2.997389078 * plaxyvr2 * plaxyvr2 + 3.020887613 * plonpar0 * plonpar0 + 
	3.044386387 * plonpar1 * plonpar1 + 3.06788516 * plonpar2 * plonpar2 + 
	3.091383696 * plonsea0 * plonsea0 + 3.114882469 * plonsea1 * plonsea1 + 
	3.138381243 * plonyvr0 * plonyvr0 + 3.161879778 * plonyvr1 * plonyvr1 + 
	3.185378551 * plonyvr7 * plonyvr7 + 3.208877325 * poakont0 * poakont0 + 
	3.23237586 * poakont1 * poakont1 + 3.255874634 * poaksea0 * poaksea0 + 
	3.279373407 * poaksea1 * poaksea1 + 3.302872181 * poaksea2 * poaksea2 + 
	3.326370716 * pontsfo0 * pontsfo0 + 3.34986949 * pontsea0 * pontsea0 + 
	3.373368263 * pontsea1 * pontsea1 + 3.396866798 * pontsea2 * pontsea2 + 
	3.420365572 * pparsea0 * pparsea0 + 3.443864346 * pparsea1 * pparsea1 + 
	3.467362881 * pparyvr0 * pparyvr0 + 3.490861654 * pparyvr1 * pparyvr1 + 
	3.514360428 * pparyvr7 * pparyvr7 + 3.537858963 * prnosea0 * prnosea0 + 
	3.561357737 * prnosea1 * prnosea1 + 3.58485651 * prnoyvr0 * prnoyvr0 + 
	3.608355045 * pseasfo0 * pseasfo0 + 3.631853819 * pseasfo1 * pseasfo1 + 
	3.655352592 * pseasfo2 * pseasfo2 + 3.678851128 * pseasfo3 * pseasfo3 + 
	3.702349901 * pseasfo4 * pseasfo4 + 3.725848675 * pseasfo5 * pseasfo5 + 
	3.74934721 * pseasfo6 * pseasfo6 + 3.772845984 * pseatpe0 * pseatpe0 + 
	3.796344757 * pseatpe1 * pseatpe1 + 3.819843292 * pseatpe2 * pseatpe2 + 
	3.843342066 * pseatpe3 * pseatpe3 + 3.866840839 * pseatyo0 * pseatyo0 + 
	3.890339375 * pseatyo1 * pseatyo1 + 3.913838148 * pseatyo2 * pseatyo2 + 
	3.937336922 * pseatyo3 * pseatyo3 + 3.960835457 * pseayvr0 * pseayvr0 + 
	3.98433423 * pseayvr1 * pseayvr1 + 4.007833004 * pseayvr2 * pseayvr2 + 
	4.031331539 * pseayvr3 * pseayvr3 + 4.054830074 * pseayvr4 * pseayvr4 + 
	4.078329086 * pseayvr5 * pseayvr5 + 4.101827621 * pseayvr6 * pseayvr6 + 
	4.125326157 * pseayvr7 * pseayvr7 + 4.148825169 * pseayvr8 * pseayvr8 + 
	4.172323704 * psfotpe0 * psfotpe0 + 4.195822239 * psfotpe1 * psfotpe1 + 
	4.219321251 * psfotpe2 * psfotpe2 + 4.242819786 * psfotpe8 * psfotpe8 + 
	4.266318321 * psfotyo0 * psfotyo0 + 4.289817333 * psfotyo1 * psfotyo1 + 
	4.313315868 * psfotyo2 * psfotyo2 + 4.336814404 * psfotyo8 * psfotyo8 + 
	4.360313416 * psfoyvr0 * psfoyvr0 + 4.383811951 * psfoyvr1 * psfoyvr1 + 
	4.407310486 * ptpetyo0 * ptpetyo0 + 4.430809498 * ptpetyo1 * ptpetyo1 + 
	4.454308033 * ptpetyo2 * ptpetyo2 + 4.477806568 * ptpetyo3 * ptpetyo3 + 
	4.50130558 * ptpeyvr0 * ptpeyvr0 + 4.524804115 * ptyoyvr0 * ptyoyvr0 + 
	4.54830265 * pyulyvr0 * pyulyvr0 + 4.571801662 * pyulyvr1 * pyulyvr1 + 
	4.595300198 * pyulyvr2 * pyulyvr2 + 4.618798733 * pyulyvr3 * pyulyvr3 + 
	4.642297745 * pyulywg0 * pyulywg0 + 4.66579628 * pyulywg1 * pyulywg1 + 
	4.689294815 * pyulywg2 * pyulywg2 + 4.712793827 * pyulywg3 * pyulywg3 + 
	4.736292362 * pyulyyz0 * pyulyyz0 + 4.759790897 * pyulyyz1 * pyulyyz1 + 
	4.783289909 * pyulyyz2 * pyulyyz2 + 4.806788445 * pyulyyz3 * pyulyyz3 + 
	4.83028698 * pyulyyz4 * pyulyyz4 + 4.853785992 * pyvrywg0 * pyvrywg0 + 
	4.877284527 * pyvrywg1 * pyvrywg1 + 4.900783062 * pyvrywg2 * pyvrywg2 + 
	4.924282074 * pyvryyz0 * pyvryyz0 + 4.947780609 * pyvryyz1 * pyvryyz1 + 
	4.971279144 * pyvryyz2 * pyvryyz2 + 4.994778156 * pywgyyz0 * pywgyyz0 + 
	5.018276691 * pywgyyz1 * pywgyyz1 + 5.041775227 * pywgyyz2 * pywgyyz2 + 
	5.065274239 * pywgyyz3 * pywgyyz3 + 5.088772774 * pbosoak0 * pbosoak0 + 
	5.112271309 * pbosoak6 * pbosoak6 + 5.135770321 * pbosbur1 * pbosbur1 + 
	5.159268856 * pbosbur2 * pbosbur2 + 5.182767391 * pbosont1 * pbosont1 + 
	5.206266403 * pbosont2 * pbosont2 + 5.229764938 * pburyvr1 * pburyvr1 + 
	5.253263474 * pburtyo1 * pburtyo1 + 5.276762486 * pburtpe1 * pburtpe1 + 
	5.300261021 * pburhnl0 * pburhnl0 + 5.323759556 * pburhnl6 * pburhnl6 + 
	5.347258568 * phnloak0 * phnloak0 + 5.370757103 * phnloak1 * phnloak1 + 
	5.394255638 * phnloak2 * phnloak2 + 5.41775465 * phnloak8 * phnloak8 + 
	5.441253185 * phnlont0 * phnlont0 + 5.46475172 * phnlont6 * phnlont6 + 
	5.488250732 * phnlywg1 * phnlywg1 + 5.511749268 * phnlyyz1 * phnlyyz1 + 
	5.53524828 * phnlyul1 * phnlyul1 + 5.558746815 * plastyo1 * plastyo1 + 
	5.58224535 * plastpe1 * plastpe1 + 5.605744362 * plaxlon0 * plaxlon0 + 
	5.629242897 * plaxlon6 * plaxlon6 + 5.652741432 * plaxlon7 * plaxlon7 + 
	5.676240444 * plaxpar0 * plaxpar0 + 5.699738979 * plaxpar6 * plaxpar6 + 
	5.723237514 * plaxpar7 * plaxpar7 + 5.746736526 * pburlon1 * pburlon1 + 
	5.770235062 * pburpar1 * pburpar1 + 5.793733597 * plonont1 * plonont1 + 
	5.817232609 * plonoak1 * plonoak1 + 5.840731144 * poakpar1 * poakpar1 + 
	5.864229679 * poaktyo1 * poaktyo1 + 5.887728691 * poaktpe1 * poaktpe1 + 
	5.911227226 * pontpar1 * pontpar1 + 5.934725761 * ponttyo1 * ponttyo1 + 
	5.958224773 * ponttpe1 * ponttpe1 + 5.981723309 * pparsfo1 * pparsfo1 + 
	6.005221844 * prnotyo1 * prnotyo1 + 6.028720856 * prnotpe1 * prnotpe1 + 
	6.052219391 * ptpeywg1 * ptpeywg1 + 6.075717926 * ptpeyyz1 * ptpeyyz1 + 
	6.099216938 * ptpeyul1 * ptpeyul1 + 6.122715473 * ptyoyul1 * ptyoyul1 + 
	6.146214008 * ptyoyyz1 * ptyoyyz1 + 6.16971302 * ptyoywg1 * ptyoywg1 + 
	6.193211555 * plaxont0 * plaxont0 + 6.216710091 * grdtimo1 * grdtimo1 + 
	6.240209103 * grdtimn1 * grdtimn1 + 6.263707638 * grdtimo2 * grdtimo2 + 
	6.287206173 * grdtimn2 * grdtimn2 + 6.310705185 * grdtimo3 * grdtimo3 + 
	6.33420372 * grdtimn3 * grdtimn3 + 6.357702255 * grdtimo4 * grdtimo4 + 
	6.381201267 * grdtimn4 * grdtimn4 + 6.404699802 * grdtimo5 * grdtimo5 + 
	6.428198338 * grdtimn5 * grdtimn5 + 6.45169735 * grdtimo6 * grdtimo6 + 
	6.475195885 * grdtimn6 * grdtimn6 + 6.49869442 * n1001ac1 * n1001ac1 + 
	6.522193432 * n1001ac2 * n1001ac2 + 6.545691967 * n1001ac3 * n1001ac3 + 
	6.569190502 * n1002ac1 * n1002ac1 + 6.592689514 * n1002ac2 * n1002ac2 + 
	6.616188049 * n1002ac3 * n1002ac3 + 6.639686584 * n1003ac1 * n1003ac1 + 
	6.663185596 * n1003ac2 * n1003ac2 + 6.686684132 * n1003ac3 * n1003ac3 + 
	6.710182667 * n1004ac1 * n1004ac1 + 6.733681679 * n1004ac2 * n1004ac2 + 
	6.757180214 * n1004ac3 * n1004ac3 + 6.780678749 * n1005ac3 * n1005ac3 + 
	6.804177761 * n1105ac3 * n1105ac3 + 6.827676296 * n1006ac3 * n1006ac3 + 
	6.851174831 * n1007ac1 * n1007ac1 + 6.874673843 * n1007ac2 * n1007ac2 + 
	6.898172379 * n1007ac3 * n1007ac3 + 6.921670914 * n1008ac1 * n1008ac1 + 
	6.945169926 * n1008ac2 * n1008ac2 + 6.968668461 * n1008ac3 * n1008ac3 + 
	6.992166996 * n1008ac4 * n1008ac4 + 7.015666008 * n1008ac5 * n1008ac5 + 
	7.039164543 * n1008ac6 * n1008ac6 + 7.062663078 * n1009ac1 * n1009ac1 + 
	7.08616209 * n1009ac2 * n1009ac2 + 7.109660625 * n1009ac3 * n1009ac3 + 
	7.133159161 * n1009ac4 * n1009ac4 + 7.156658173 * n1009ac5 * n1009ac5 + 
	7.180156708 * n1010ac1 * n1010ac1 + 7.203655243 * n1010ac2 * n1010ac2 + 
	7.227154255 * n1010ac3 * n1010ac3 + 7.25065279 * n1010ac4 * n1010ac4 + 
	7.274151325 * n1010ac5 * n1010ac5 + 7.297650337 * n1010ac6 * n1010ac6 + 
	7.321148872 * n1011ac1 * n1011ac1 + 7.344647408 * n1011ac2 * n1011ac2 + 
	7.36814642 * n1011ac3 * n1011ac3 + 7.391644955 * n1011ac4 * n1011ac4 + 
	7.41514349 * n1011ac5 * n1011ac5 + 7.438642502 * n1011ac6 * n1011ac6 + 
	7.462141037 * n1012ac1 * n1012ac1 + 7.485639572 * n1012ac2 * n1012ac2 + 
	7.509138584 * n1012ac3 * n1012ac3 + 7.532637119 * n1012ac4 * n1012ac4 + 
	7.556135654 * n1012ac5 * n1012ac5 + 7.579634666 * n1013ac3 * n1013ac3 + 
	7.603133202 * n1013ac4 * n1013ac4 + 7.626631737 * n1013ac5 * n1013ac5 + 
	7.650130749 * n1013ac6 * n1013ac6 + 7.673629284 * n1014ac3 * n1014ac3 + 
	7.697127819 * n1014ac4 * n1014ac4 + 7.720626831 * n1014ac5 * n1014ac5 + 
	7.744125366 * n1014ac6 * n1014ac6 + 7.767623901 * n1015ac3 * n1015ac3 + 
	7.791122913 * n1015ac4 * n1015ac4 + 7.814621449 * n1015ac5 * n1015ac5 + 
	7.838119984 * n1015ac6 * n1015ac6 + 7.861618996 * n1016ac3 * n1016ac3 + 
	7.885117531 * n1016ac4 * n1016ac4 + 7.9086160660000004 * n1016ac5 * n1016ac5 + 
	7.932115078 * n1016ac6 * n1016ac6 + 7.955613613 * n1017ac3 * n1017ac3 + 
	7.979112148 * n1017ac4 * n1017ac4 + 8.00261116 * n1017ac5 * n1017ac5 + 
	8.026109695 * n1017ac6 * n1017ac6 + 8.049608231 * n1018ac1 * n1018ac1 + 
	8.073106766 * n1018ac2 * n1018ac2 + 8.096605301 * n1018ac3 * n1018ac3 + 
	8.12010479 * n1018ac4 * n1018ac4 + 8.143603325 * n1018ac5 * n1018ac5 + 
	8.16710186 * n1018ac6 * n1018ac6 + 8.190600395 * n1019ac1 * n1019ac1 + 
	8.21409893 * n1019ac2 * n1019ac2 + 8.237597466 * n1019ac3 * n1019ac3 + 
	8.261096954 * n1019ac4 * n1019ac4 + 8.28459549 * n1019ac5 * n1019ac5 + 
	8.308094025 * n1020ac1 * n1020ac1 + 8.33159256 * n1020ac2 * n1020ac2 + 
	8.355091095 * n1020ac3 * n1020ac3 + 8.37858963 * n1020ac4 * n1020ac4 + 
	8.402089119 * n1020ac5 * n1020ac5 + 8.425587654 * n1020ac6 * n1020ac6 + 
	8.449086189 * n1021ac1 * n1021ac1 + 8.472584724 * n1021ac2 * n1021ac2 + 
	8.49608326 * n1021ac3 * n1021ac3 + 8.519581795 * n1021ac4 * n1021ac4 + 
	8.543081284 * n1021ac5 * n1021ac5 + 8.566579819 * n1022ac1 * n1022ac1 + 
	8.590078354 * n1023ac1 * n1023ac1 + 8.613576889 * n1026ac1 * n1026ac1 + 
	8.637075424 * n1027ac1 * n1027ac1 + 8.660573959 * n1028ac1 * n1028ac1 + 
	8.684073448 * n1029ac1 * n1029ac1 + 8.707571983 * n1030ac1 * n1030ac1 + 
	8.731070518 * n1032ac1 * n1032ac1 + 8.754569054 * n1032ac2 * n1032ac2 + 
	8.778067589 * n1032ac3 * n1032ac3 + 8.801566124 * n1032ac4 * n1032ac4 + 
	8.825065613 * n1032ac5 * n1032ac5 + 8.848564148 * n1033ac1 * n1033ac1 + 
	8.872062683 * n1033ac2 * n1033ac2 + 8.895561218 * n1033ac3 * n1033ac3 + 
	8.919059753 * n1033ac4 * n1033ac4 + 8.942558289 * n1033ac5 * n1033ac5 + 
	8.966057777 * n1034ac1 * n1034ac1 + 8.989556313 * n1034ac2 * n1034ac2 + 
	9.013054848 * n1034ac3 * n1034ac3 + 9.036553383 * n1035ac1 * n1035ac1 + 
	9.060051918 * n1035ac2 * n1035ac2 + 9.083550453 * n1035ac3 * n1035ac3 + 
	9.107049942 * n1035ac4 * n1035ac4 + 9.130548477 * n1035ac5 * n1035ac5 + 
	9.154047012 * n1036ac1 * n1036ac1 + 9.177545547 * n1036ac2 * n1036ac2 + 
	9.201044083 * n1036ac3 * n1036ac3 + 9.224542618 * n1037ac4 * n1037ac4 + 
	9.248042107 * n1037ac5 * n1037ac5 + 9.271540642 * n1038ac4 * n1038ac4 + 
	9.295039177 * n1038ac5 * n1038ac5 + 9.318537712 * n1039ac4 * n1039ac4 + 
	9.342036247 * n1039ac5 * n1039ac5 + 9.365534782 * n1040ac4 * n1040ac4 + 
	9.389034271 * n1040ac5 * n1040ac5 + 9.412532806 * n1040ac6 * n1040ac6 + 
	9.436031342 * n1041ac4 * n1041ac4 + 9.459529877 * n1041ac5 * n1041ac5 + 
	9.483028412 * n1041ac6 * n1041ac6 + 9.506526947 * n1042ac4 * n1042ac4 + 
	9.530026436 * n1042ac5 * n1042ac5 + 9.553524971 * n1042ac6 * n1042ac6 + 
	9.577023506 * n1043ac1 * n1043ac1 + 9.600522041 * n1043ac2 * n1043ac2 + 
	9.624020576 * n1043ac3 * n1043ac3 + 9.647519112 * n1044ac1 * n1044ac1 + 
	9.6710186 * n1044ac2 * n1044ac2 + 9.694517136 * n1044ac3 * n1044ac3 + 
	9.718015671 * n1046ac3 * n1046ac3 + 9.741514206 * n1047ac1 * n1047ac1 + 
	9.765012741 * n1047ac2 * n1047ac2 + 9.788511276 * n1047ac3 * n1047ac3 + 
	9.812010765 * n1050ac3 * n1050ac3 + 9.8355093 * n1050ac4 * n1050ac4 + 
	9.859007835 * n1050ac5 * n1050ac5 + 9.882506371 * n1051ac1 * n1051ac1 + 
	9.906004906 * n1051ac2 * n1051ac2 + 9.929503441 * n1051ac3 * n1051ac3 + 
	9.95300293 * n1051ac4 * n1051ac4 + 9.976501465 * n1051ac5 * n1051ac5 + 10.0 * 
	n1051ac6 * n1051ac6 - 0.23519*pboshnl0 - 0.23519*pboshnl1 + 0.02469*pboshnl7 + 
	0.02665*pboshnl8 - 0.13092*pboslax0 - 0.13092*pboslax1 + 0.03462*pboslax7 - 
	0.12596*pbossea0 - 0.12596*pbossea1 - 0.12596*pbossea2 - 0.13481*pbossfo0 - 
	0.13481*pbossfo1 + 0.05629*pbostpe1 + 0.0528*pbostpe2 + 0.04525*pbostyo1 + 
	0.04206*pbostyo2 - 0.02799*pbosyul0 - 0.02799*pbosyul1 - 0.02799*pbosyul2 - 
	0.02799*pbosyul3 - 0.02799*pbosyul4 - 0.12674*pbosyvr0 - 0.12674*pbosyvr1 - 
	0.12674*pbosyvr2 - 0.07408*pbosywg0 - 0.07408*pbosywg1 + 0.01789*pbosywg7 - 
	0.03418*pbosyyz0 - 0.03418*pbosyyz1 - 0.03032*pburoak0 - 0.03032*pburoak1 - 
	0.05527*pbursea0 - 0.05527*pbursea1 - 0.03033*pbursfo0 - 0.12871*phnllax0 - 
	0.12871*phnllax1 - 0.12871*phnllax2 - 0.12871*phnllax3 - 0.32504*phnllon0 + 
	0.03143*phnllon6 - 0.3343*phnlpar0 + 0.03166*phnlpar6 - 0.13391*phnlsea0 - 
	0.13391*phnlsea1 - 0.13391*phnlsea2 - 0.12174*phnlsfo0 - 0.12174*phnlsfo1 + 
	0.0377*phnlsfo7 - 0.1351*phnlyvr0 - 0.1351*phnlyvr1 + 0.01809*phnlyvr7 - 
	0.05205*plassea0 - 0.05205*plassea1 - 0.05774*plasyvr0 + 0.01359*plasyvr6 - 
	0.03072*plaxoak0 - 0.03072*plaxoak1 - 0.03072*plaxoak2 - 0.05607*plaxsea0 - 
	0.05607*plaxsea1 - 0.05607*plaxsea2 - 0.05607*plaxsea3 - 0.05607*plaxsea4 - 
	0.05607*plaxsea5 - 0.05607*plaxsea6 - 0.05607*plaxsea7 - 0.05607*plaxsea8 - 
	0.05607*plaxsea9 - 0.03072*plaxsfo0 - 0.03072*plaxsfo1 - 0.03072*plaxsfo2 - 
	0.03072*plaxsfo3 - 0.03072*plaxsfo4 - 0.03072*plaxsfo5 - 0.03072*plaxsfo6 - 
	0.30609*plaxtpe0 - 0.30609*plaxtpe1 - 0.30609*plaxtpe2 + 0.02571*plaxtpe8 - 
	0.2516*plaxtyo0 - 0.2516*plaxtyo1 - 0.2516*plaxtyo2 + 0.0273*plaxtyo8 - 
	0.06181*plaxyvr0 - 0.06181*plaxyvr1 - 0.06181*plaxyvr2 - 0.02622*plonpar0 - 
	0.02622*plonpar1 - 0.02622*plonpar2 - 0.22255*plonsea0 - 0.22255*plonsea1 - 
	0.21935*plonyvr0 - 0.21935*plonyvr1 + 0.02248*plonyvr7 - 0.03148*poakont0 - 
	0.03148*poakont1 - 0.04311*poaksea0 - 0.04311*poaksea1 - 0.04311*poaksea2 - 
	0.03152*pontsfo0 - 0.05614*pontsea0 - 0.05614*pontsea1 - 0.05614*pontsea2 - 
	0.23205*pparsea0 - 0.23205*pparsea1 - 0.22885*pparyvr0 - 0.22885*pparyvr1 + 
	0.02247*pparyvr7 - 0.0382*prnosea0 - 0.0382*prnosea1 - 0.04394*prnoyvr0 - 
	0.04347*pseasfo0 - 0.04347*pseasfo1 - 0.04347*pseasfo2 - 0.04347*pseasfo3 - 
	0.04347*pseasfo4 - 0.04347*pseasfo5 - 0.04347*pseasfo6 - 0.27573*pseatpe0 - 
	0.27573*pseatpe1 - 0.27573*pseatpe2 - 0.27573*pseatpe3 - 0.22283*pseatyo0 - 
	0.22283*pseatyo1 - 0.22283*pseatyo2 - 0.22283*pseatyo3 - 0.01928*pseayvr0 - 
	0.01928*pseayvr1 - 0.01928*pseayvr2 - 0.01928*pseayvr3 - 0.01928*pseayvr4 - 
	0.01928*pseayvr5 - 0.01928*pseayvr6 - 0.01928*pseayvr7 - 0.01928*pseayvr8 - 
	0.29218*psfotpe0 - 0.29218*psfotpe1 - 0.29218*psfotpe2 + 0.02702*psfotpe8 - 
	0.2378*psfotyo0 - 0.2378*psfotyo1 - 0.2378*psfotyo2 + 0.0285*psfotyo8 - 
	0.04898*psfoyvr0 - 0.04898*psfoyvr1 - 0.07188*ptpetyo0 - 0.07188*ptpetyo1 - 
	0.07188*ptpetyo2 - 0.07188*ptpetyo3 - 0.27146*ptpeyvr0 - 0.21887*ptyoyvr0 - 
	0.11677*pyulyvr0 - 0.11677*pyulyvr1 - 0.11677*pyulyvr2 - 0.11677*pyulyvr3 - 
	0.06398*pyulywg0 - 0.06398*pyulywg1 - 0.06398*pyulywg2 - 0.06398*pyulywg3 - 
	0.02995*pyulyyz0 - 0.02995*pyulyyz1 - 0.02995*pyulyyz2 - 0.02995*pyulyyz3 - 
	0.02995*pyulyyz4 - 0.06529*pyvrywg0 - 0.06529*pyvrywg1 - 0.06529*pyvrywg2 - 
	0.10723*pyvryyz0 - 0.10723*pyvryyz1 - 0.10723*pyvryyz2 - 0.05507*pywgyyz0 - 
	0.05507*pywgyyz1 - 0.05507*pywgyyz2 - 0.05507*pywgyyz3 - 0.13432*pbosoak0 + 
	0.01157*pbosoak6 + 0.03469*pbosbur1 + 0.04664*pbosbur2 + 0.03725*pbosont1 + 
	0.04996*pbosont2 + 0.01353*pburyvr1 + 0.02686*pburtyo1 + 0.02529*pburtpe1 - 
	0.12896*pburhnl0 + 0.01132*pburhnl6 - 0.12225*phnloak0 - 0.12225*phnloak1 - 
	0.12225*phnloak2 + 0.03719*phnloak8 - 0.13063*phnlont0 + 0.0116*phnlont6 + 
	0.01893*phnlywg1 + 0.02585*phnlyyz1 + 0.02429*phnlyul1 + 0.02076*plastyo1 + 
	0.01956*plastpe1 - 0.25019*plaxlon0 + 0.02843*plaxlon6 + 0.03096*plaxlon7 - 
	0.25943*plaxpar0 + 0.02868*plaxpar6 + 0.03123*plaxpar7 + 0.02834*pburlon1 + 
	0.02859*pburpar1 + 0.02984*plonont1 + 0.0197*plonoak1 + 0.0198*poakpar1 + 
	0.028*poaktyo1 + 0.02654*poaktpe1 + 0.03012*pontpar1 + 0.02605*ponttyo1 + 
	0.02451*ponttpe1 + 0.01968*pparsfo1 + 0.02132*prnotyo1 + 0.0202*prnotpe1 + 
	0.03219*ptpeywg1 + 0.04208*ptpeyyz1 + 0.05292*ptpeyul1 + 0.0427*ptyoyul1 + 
	0.03452*ptyoyyz1 + 0.02803*ptyoywg1 - 0.01351*plaxont0 + 0.457*grdtimo1 - 
	0.13333*grdtimn1 + 0.318*grdtimo2 - 0.10692*grdtimn2 + 0.206*grdtimo3 + 
	0.25*grdtimo4 + 0.163*grdtimo5 + 0.095*grdtimo6 + 14.0062*n1001ac1 + 
	10.44277*n1001ac2 + 7.65023*n1001ac3 + 16.89049*n1002ac1 + 12.39107*n1002ac2 + 
	8.97425*n1002ac3 + 13.46672*n1003ac1 + 10.03075*n1003ac2 + 7.34337*n1003ac3 + 
	17.27667*n1004ac1 + 12.68603*n1004ac2 + 9.19391*n1004ac3 + 9.1777*n1005ac3 + 
	9.19396*n1105ac3 + 9.3651*n1006ac3 + 27.19154*n1007ac1 + 20.25861*n1007ac2 + 
	14.83356*n1007ac3 + 11.99082*n1008ac1 + 8.39427*n1008ac2 + 5.87068*n1008ac3 + 
	3.99147*n1008ac4 + 4.41133*n1008ac5 + 3.43737*n1008ac6 + 9.38535*n1009ac1 + 
	6.65893*n1009ac2 + 4.70526*n1009ac3 + 3.27055*n1009ac4 + 3.58556*n1009ac5 + 
	2.8843*n1010ac1 + 1.9483*n1010ac2 + 1.32401*n1010ac3 + 0.84307*n1010ac4 + 
	0.95496*n1010ac5 + 0.67624*n1010ac6 + 9.10651*n1011ac1 + 6.44597*n1011ac2 + 
	4.54666*n1011ac3 + 3.1484*n1011ac4 + 3.45638*n1011ac5 + 2.76114*n1011ac6 + 
	6.50105*n1012ac1 + 4.71063*n1012ac2 + 3.38125*n1012ac3 + 2.42747*n1012ac4 + 
	2.6306*n1012ac5 + 4.58479*n1013ac3 + 3.17776*n1013ac4 + 3.48744*n1013ac5 + 
	2.78943*n1013ac6 + 4.4958*n1014ac3 + 3.10922*n1014ac4 + 3.41495*n1014ac5 + 
	2.7234*n1014ac6 + 4.60806*n1015ac3 + 3.19568*n1015ac4 + 3.50639*n1015ac5 + 
	2.80669*n1015ac6 + 4.51701*n1016ac3 + 3.12556*n1016ac4 + 3.43223*n1016ac5 + 
	2.73914*n1016ac6 + 4.52672*n1017ac3 + 3.13304*n1017ac4 + 3.44014*n1017ac5 + 
	2.74634*n1017ac6 + 7.67623*n1018ac1 + 5.35357*n1018ac2 + 3.7331*n1018ac3 + 
	2.52182*n1018ac4 + 2.79372*n1018ac5 + 2.15753*n1018ac6 + 9.00153*n1019ac1 + 
	6.36578*n1019ac2 + 4.48694*n1019ac3 + 3.1024*n1019ac4 + 3.40774*n1019ac5 + 
	4.79193*n1020ac1 + 3.40527*n1020ac2 + 2.40909*n1020ac3 + 1.67875*n1020ac4 + 
	1.83876*n1020ac5 + 1.48129*n1020ac6 + 6.11723*n1021ac1 + 4.41748*n1021ac2 + 
	3.16293*n1021ac3 + 2.25933*n1021ac4 + 2.45278*n1021ac5 + 40.32039*n1022ac1 + 
	29.4082*n1023ac1 + 33.02493*n1026ac1 + 37.74599*n1027ac1 + 36.5415*n1028ac1 + 
	42.82423*n1029ac1 + 40.35146*n1030ac1 + 23.8204*n1032ac1 + 16.92001*n1032ac2 + 
	11.96625*n1032ac3 + 8.33277*n1032ac4 + 9.12932*n1032ac5 + 18.08295*n1033ac1 + 
	13.0472*n1033ac2 + 9.33593*n1033ac3 + 6.66027*n1033ac4 + 7.23383*n1033ac5 + 
	15.09963*n1034ac1 + 11.02328*n1034ac2 + 7.95559*n1034ac3 + 13.78819*n1035ac1 + 
	10.02165*n1035ac2 + 7.20963*n1035ac3 + 5.19931*n1035ac4 + 5.62541*n1035ac5 + 
	15.74868*n1036ac1 + 11.51899*n1036ac2 + 8.32477*n1036ac3 + 3.5884*n1037ac4 + 
	3.92172*n1037ac5 + 2.38548*n1038ac4 + 2.5862*n1038ac5 + 2.75896*n1039ac4 + 
	2.98119*n1039ac5 + 2.29039*n1040ac4 + 2.54896*n1040ac5 + 1.93457*n1040ac6 + 
	1.20291*n1041ac4 + 1.33552*n1041ac5 + 1.02289*n1041ac6 + 1.08748*n1042ac4 + 
	1.21344*n1042ac5 + 0.91168*n1042ac6 + 13.18536*n1043ac1 + 9.81585*n1043ac2 + 
	7.18333*n1043ac3 + 14.09924*n1044ac1 + 10.51384*n1044ac2 + 7.70315*n1044ac3 + 
	9.52503*n1046ac3 + 31.37592*n1047ac1 + 23.19986*n1047ac2 + 16.89706*n1047ac3 + 
	4.50544*n1050ac3 + 3.11665*n1050ac4 + 3.42281*n1050ac5 + 3.80996*n1051ac1 + 
	2.65528*n1051ac2 + 1.85054*n1051ac3 + 1.24858*n1051ac4 + 1.38382*n1051ac5 + 
	1.06688*n1051ac6;

subject to revenues:
	0 <= 0.23519*pboshnl0 + 0.23519*pboshnl1 - 0.02469*pboshnl7 - 0.02665*pboshnl8 
	+ 0.13092*pboslax0 + 0.13092*pboslax1 - 0.03462*pboslax7 + 0.12596*pbossea0 + 
	0.12596*pbossea1 + 0.12596*pbossea2 + 0.13481*pbossfo0 + 0.13481*pbossfo1 - 
	0.05629*pbostpe1 - 0.0528*pbostpe2 - 0.04525*pbostyo1 - 0.04206*pbostyo2 + 
	0.02799*pbosyul0 + 0.02799*pbosyul1 + 0.02799*pbosyul2 + 0.02799*pbosyul3 + 
	0.02799*pbosyul4 + 0.12674*pbosyvr0 + 0.12674*pbosyvr1 + 0.12674*pbosyvr2 + 
	0.07408*pbosywg0 + 0.07408*pbosywg1 - 0.01789*pbosywg7 + 0.03418*pbosyyz0 + 
	0.03418*pbosyyz1 + 0.03032*pburoak0 + 0.03032*pburoak1 + 0.05527*pbursea0 + 
	0.05527*pbursea1 + 0.03033*pbursfo0 + 0.12871*phnllax0 + 0.12871*phnllax1 + 
	0.12871*phnllax2 + 0.12871*phnllax3 + 0.32504*phnllon0 - 0.03143*phnllon6 + 
	0.3343*phnlpar0 - 0.03166*phnlpar6 + 0.13391*phnlsea0 + 0.13391*phnlsea1 + 
	0.13391*phnlsea2 + 0.12174*phnlsfo0 + 0.12174*phnlsfo1 - 0.0377*phnlsfo7 + 
	0.1351*phnlyvr0 + 0.1351*phnlyvr1 - 0.01809*phnlyvr7 + 0.05205*plassea0 + 
	0.05205*plassea1 + 0.05774*plasyvr0 - 0.01359*plasyvr6 + 0.03072*plaxoak0 + 
	0.03072*plaxoak1 + 0.03072*plaxoak2 + 0.05607*plaxsea0 + 0.05607*plaxsea1 + 
	0.05607*plaxsea2 + 0.05607*plaxsea3 + 0.05607*plaxsea4 + 0.05607*plaxsea5 + 
	0.05607*plaxsea6 + 0.05607*plaxsea7 + 0.05607*plaxsea8 + 0.05607*plaxsea9 + 
	0.03072*plaxsfo0 + 0.03072*plaxsfo1 + 0.03072*plaxsfo2 + 0.03072*plaxsfo3 + 
	0.03072*plaxsfo4 + 0.03072*plaxsfo5 + 0.03072*plaxsfo6 + 0.30609*plaxtpe0 + 
	0.30609*plaxtpe1 + 0.30609*plaxtpe2 - 0.02571*plaxtpe8 + 0.2516*plaxtyo0 + 
	0.2516*plaxtyo1 + 0.2516*plaxtyo2 - 0.0273*plaxtyo8 + 0.06181*plaxyvr0 + 
	0.06181*plaxyvr1 + 0.06181*plaxyvr2 + 0.02622*plonpar0 + 0.02622*plonpar1 + 
	0.02622*plonpar2 + 0.22255*plonsea0 + 0.22255*plonsea1 + 0.21935*plonyvr0 + 
	0.21935*plonyvr1 - 0.02248*plonyvr7 + 0.03148*poakont0 + 0.03148*poakont1 + 
	0.04311*poaksea0 + 0.04311*poaksea1 + 0.04311*poaksea2 + 0.03152*pontsfo0 + 
	0.05614*pontsea0 + 0.05614*pontsea1 + 0.05614*pontsea2 + 0.23205*pparsea0 + 
	0.23205*pparsea1 + 0.22885*pparyvr0 + 0.22885*pparyvr1 - 0.02247*pparyvr7 + 
	0.0382*prnosea0 + 0.0382*prnosea1 + 0.04394*prnoyvr0 + 0.04347*pseasfo0 + 
	0.04347*pseasfo1 + 0.04347*pseasfo2 + 0.04347*pseasfo3 + 0.04347*pseasfo4 + 
	0.04347*pseasfo5 + 0.04347*pseasfo6 + 0.27573*pseatpe0 + 0.27573*pseatpe1 + 
	0.27573*pseatpe2 + 0.27573*pseatpe3 + 0.22283*pseatyo0 + 0.22283*pseatyo1 + 
	0.22283*pseatyo2 + 0.22283*pseatyo3 + 0.01928*pseayvr0 + 0.01928*pseayvr1 + 
	0.01928*pseayvr2 + 0.01928*pseayvr3 + 0.01928*pseayvr4 + 0.01928*pseayvr5 + 
	0.01928*pseayvr6 + 0.01928*pseayvr7 + 0.01928*pseayvr8 + 0.29218*psfotpe0 + 
	0.29218*psfotpe1 + 0.29218*psfotpe2 - 0.02702*psfotpe8 + 0.2378*psfotyo0 + 
	0.2378*psfotyo1 + 0.2378*psfotyo2 - 0.0285*psfotyo8 + 0.04898*psfoyvr0 + 
	0.04898*psfoyvr1 + 0.07188*ptpetyo0 + 0.07188*ptpetyo1 + 0.07188*ptpetyo2 + 
	0.07188*ptpetyo3 + 0.27146*ptpeyvr0 + 0.21887*ptyoyvr0 + 0.11677*pyulyvr0 + 
	0.11677*pyulyvr1 + 0.11677*pyulyvr2 + 0.11677*pyulyvr3 + 0.06398*pyulywg0 + 
	0.06398*pyulywg1 + 0.06398*pyulywg2 + 0.06398*pyulywg3 + 0.02995*pyulyyz0 + 
	0.02995*pyulyyz1 + 0.02995*pyulyyz2 + 0.02995*pyulyyz3 + 0.02995*pyulyyz4 + 
	0.06529*pyvrywg0 + 0.06529*pyvrywg1 + 0.06529*pyvrywg2 + 0.10723*pyvryyz0 + 
	0.10723*pyvryyz1 + 0.10723*pyvryyz2 + 0.05507*pywgyyz0 + 0.05507*pywgyyz1 + 
	0.05507*pywgyyz2 + 0.05507*pywgyyz3 + 0.13432*pbosoak0 - 0.01157*pbosoak6 - 
	0.03469*pbosbur1 - 0.04664*pbosbur2 - 0.03725*pbosont1 - 0.04996*pbosont2 - 
	0.01353*pburyvr1 - 0.02686*pburtyo1 - 0.02529*pburtpe1 + 0.12896*pburhnl0 - 
	0.01132*pburhnl6 + 0.12225*phnloak0 + 0.12225*phnloak1 + 0.12225*phnloak2 - 
	0.03719*phnloak8 + 0.13063*phnlont0 - 0.0116*phnlont6 - 0.01893*phnlywg1 - 
	0.02585*phnlyyz1 - 0.02429*phnlyul1 - 0.02076*plastyo1 - 0.01956*plastpe1 + 
	0.25019*plaxlon0 - 0.02843*plaxlon6 - 0.03096*plaxlon7 + 0.25943*plaxpar0 - 
	0.02868*plaxpar6 - 0.03123*plaxpar7 - 0.02834*pburlon1 - 0.02859*pburpar1 - 
	0.02984*plonont1 - 0.0197*plonoak1 - 0.0198*poakpar1 - 0.028*poaktyo1 - 
	0.02654*poaktpe1 - 0.03012*pontpar1 - 0.02605*ponttyo1 - 0.02451*ponttpe1 - 
	0.01968*pparsfo1 - 0.02132*prnotyo1 - 0.0202*prnotpe1 - 0.03219*ptpeywg1 - 
	0.04208*ptpeyyz1 - 0.05292*ptpeyul1 - 0.0427*ptyoyul1 - 0.03452*ptyoyyz1 - 
	0.02803*ptyoywg1 + 0.01351*plaxont0;
subject to acocosts:
	0 <= 0.457*grdtimo1 - 0.13333*grdtimn1 + 0.318*grdtimo2 - 0.10692*grdtimn2 + 
	0.206*grdtimo3 + 0.25*grdtimo4 + 0.163*grdtimo5 + 0.095*grdtimo6 + 
	14.0062*n1001ac1 + 10.44277*n1001ac2 + 7.65023*n1001ac3 + 16.89049*n1002ac1 + 
	12.39107*n1002ac2 + 8.97425*n1002ac3 + 13.46672*n1003ac1 + 10.03075*n1003ac2 + 
	7.34337*n1003ac3 + 17.27667*n1004ac1 + 12.68603*n1004ac2 + 9.19391*n1004ac3 + 
	9.1777*n1005ac3 + 9.19396*n1105ac3 + 9.3651*n1006ac3 + 27.19154*n1007ac1 + 
	20.25861*n1007ac2 + 14.83356*n1007ac3 + 11.99082*n1008ac1 + 8.39427*n1008ac2 + 
	5.87068*n1008ac3 + 3.99147*n1008ac4 + 4.41133*n1008ac5 + 3.43737*n1008ac6 + 
	9.38535*n1009ac1 + 6.65893*n1009ac2 + 4.70526*n1009ac3 + 3.27055*n1009ac4 + 
	3.58556*n1009ac5 + 2.8843*n1010ac1 + 1.9483*n1010ac2 + 1.32401*n1010ac3 + 
	0.84307*n1010ac4 + 0.95496*n1010ac5 + 0.67624*n1010ac6 + 9.10651*n1011ac1 + 
	6.44597*n1011ac2 + 4.54666*n1011ac3 + 3.1484*n1011ac4 + 3.45638*n1011ac5 + 
	2.76114*n1011ac6 + 6.50105*n1012ac1 + 4.71063*n1012ac2 + 3.38125*n1012ac3 + 
	2.42747*n1012ac4 + 2.6306*n1012ac5 + 4.58479*n1013ac3 + 3.17776*n1013ac4 + 
	3.48744*n1013ac5 + 2.78943*n1013ac6 + 4.4958*n1014ac3 + 3.10922*n1014ac4 + 
	3.41495*n1014ac5 + 2.7234*n1014ac6 + 4.60806*n1015ac3 + 3.19568*n1015ac4 + 
	3.50639*n1015ac5 + 2.80669*n1015ac6 + 4.51701*n1016ac3 + 3.12556*n1016ac4 + 
	3.43223*n1016ac5 + 2.73914*n1016ac6 + 4.52672*n1017ac3 + 3.13304*n1017ac4 + 
	3.44014*n1017ac5 + 2.74634*n1017ac6 + 7.67623*n1018ac1 + 5.35357*n1018ac2 + 
	3.7331*n1018ac3 + 2.52182*n1018ac4 + 2.79372*n1018ac5 + 2.15753*n1018ac6 + 
	9.00153*n1019ac1 + 6.36578*n1019ac2 + 4.48694*n1019ac3 + 3.1024*n1019ac4 + 
	3.40774*n1019ac5 + 4.79193*n1020ac1 + 3.40527*n1020ac2 + 2.40909*n1020ac3 + 
	1.67875*n1020ac4 + 1.83876*n1020ac5 + 1.48129*n1020ac6 + 6.11723*n1021ac1 + 
	4.41748*n1021ac2 + 3.16293*n1021ac3 + 2.25933*n1021ac4 + 2.45278*n1021ac5 + 
	40.32039*n1022ac1 + 29.4082*n1023ac1 + 33.02493*n1026ac1 + 37.74599*n1027ac1 + 
	36.5415*n1028ac1 + 42.82423*n1029ac1 + 40.35146*n1030ac1 + 23.8204*n1032ac1 + 
	16.92001*n1032ac2 + 11.96625*n1032ac3 + 8.33277*n1032ac4 + 9.12932*n1032ac5 + 
	18.08295*n1033ac1 + 13.0472*n1033ac2 + 9.33593*n1033ac3 + 6.66027*n1033ac4 + 
	7.23383*n1033ac5 + 15.09963*n1034ac1 + 11.02328*n1034ac2 + 7.95559*n1034ac3 + 
	13.78819*n1035ac1 + 10.02165*n1035ac2 + 7.20963*n1035ac3 + 5.19931*n1035ac4 + 
	5.62541*n1035ac5 + 15.74868*n1036ac1 + 11.51899*n1036ac2 + 8.32477*n1036ac3 + 
	3.5884*n1037ac4 + 3.92172*n1037ac5 + 2.38548*n1038ac4 + 2.5862*n1038ac5 + 
	2.75896*n1039ac4 + 2.98119*n1039ac5 + 2.29039*n1040ac4 + 2.54896*n1040ac5 + 
	1.93457*n1040ac6 + 1.20291*n1041ac4 + 1.33552*n1041ac5 + 1.02289*n1041ac6 + 
	1.08748*n1042ac4 + 1.21344*n1042ac5 + 0.91168*n1042ac6 + 13.18536*n1043ac1 + 
	9.81585*n1043ac2 + 7.18333*n1043ac3 + 14.09924*n1044ac1 + 10.51384*n1044ac2 + 
	7.70315*n1044ac3 + 9.52503*n1046ac3 + 31.37592*n1047ac1 + 23.19986*n1047ac2 + 
	16.89706*n1047ac3 + 4.50544*n1050ac3 + 3.11665*n1050ac4 + 3.42281*n1050ac5 + 
	3.80996*n1051ac1 + 2.65528*n1051ac2 + 1.85054*n1051ac3 + 1.24858*n1051ac4 + 
	1.38382*n1051ac5 + 1.06688*n1051ac6;
subject to systdept:
	0 <= n1001ac1 + n1001ac2 + n1001ac3 + 2.0*n1002ac1 + 2.0*n1002ac2 + 
	2.0*n1002ac3 + n1003ac1 + n1003ac2 + n1003ac3 + 2.0*n1004ac1 + 2.0*n1004ac2 + 
	2.0*n1004ac3 + 2.0*n1005ac3 + 2.0*n1105ac3 + 2.0*n1006ac3 + 2.0*n1007ac1 + 
	2.0*n1007ac2 + 2.0*n1007ac3 + 3.0*n1008ac1 + 3.0*n1008ac2 + 3.0*n1008ac3 + 
	3.0*n1008ac4 + 3.0*n1008ac5 + 3.0*n1008ac6 + 2.0*n1009ac1 + 2.0*n1009ac2 + 
	2.0*n1009ac3 + 2.0*n1009ac4 + 2.0*n1009ac5 + n1010ac1 + n1010ac2 + n1010ac3 + 
	n1010ac4 + n1010ac5 + n1010ac6 + 2.0*n1011ac1 + 2.0*n1011ac2 + 2.0*n1011ac3 + 
	2.0*n1011ac4 + 2.0*n1011ac5 + 2.0*n1011ac6 + n1012ac1 + n1012ac2 + n1012ac3 + 
	n1012ac4 + n1012ac5 + 2.0*n1013ac3 + 2.0*n1013ac4 + 2.0*n1013ac5 + 2.0*n1013ac6 
	+ 2.0*n1014ac3 + 2.0*n1014ac4 + 2.0*n1014ac5 + 2.0*n1014ac6 + 2.0*n1015ac3 + 
	2.0*n1015ac4 + 2.0*n1015ac5 + 2.0*n1015ac6 + 2.0*n1016ac3 + 2.0*n1016ac4 + 
	2.0*n1016ac5 + 2.0*n1016ac6 + 2.0*n1017ac3 + 2.0*n1017ac4 + 2.0*n1017ac5 + 
	2.0*n1017ac6 + 2.0*n1018ac1 + 2.0*n1018ac2 + 2.0*n1018ac3 + 2.0*n1018ac4 + 
	2.0*n1018ac5 + 2.0*n1018ac6 + 2.0*n1019ac1 + 2.0*n1019ac2 + 2.0*n1019ac3 + 
	2.0*n1019ac4 + 2.0*n1019ac5 + n1020ac1 + n1020ac2 + n1020ac3 + n1020ac4 + 
	n1020ac5 + n1020ac6 + n1021ac1 + n1021ac2 + n1021ac3 + n1021ac4 + n1021ac5 + 
	3.0*n1022ac1 + 3.0*n1023ac1 + 3.0*n1026ac1 + 3.0*n1027ac1 + 3.0*n1028ac1 + 
	5.0*n1029ac1 + 4.0*n1030ac1 + 5.0*n1032ac1 + 5.0*n1032ac2 + 5.0*n1032ac3 + 
	5.0*n1032ac4 + 5.0*n1032ac5 + 3.0*n1033ac1 + 3.0*n1033ac2 + 3.0*n1033ac3 + 
	3.0*n1033ac4 + 3.0*n1033ac5 + 2.0*n1034ac1 + 2.0*n1034ac2 + 2.0*n1034ac3 + 
	2.0*n1035ac1 + 2.0*n1035ac2 + 2.0*n1035ac3 + 2.0*n1035ac4 + 2.0*n1035ac5 + 
	2.0*n1036ac1 + 2.0*n1036ac2 + 2.0*n1036ac3 + 2.0*n1037ac4 + 2.0*n1037ac5 + 
	n1038ac4 + n1038ac5 + n1039ac4 + n1039ac5 + 2.0*n1040ac4 + 2.0*n1040ac5 + 
	2.0*n1040ac6 + n1041ac4 + n1041ac5 + n1041ac6 + n1042ac4 + n1042ac5 + n1042ac6 
	+ n1043ac1 + n1043ac2 + n1043ac3 + n1044ac1 + n1044ac2 + n1044ac3 + 
	2.0*n1046ac3 + 3.0*n1047ac1 + 3.0*n1047ac2 + 3.0*n1047ac3 + 2.0*n1050ac3 + 
	2.0*n1050ac4 + 2.0*n1050ac5 + n1051ac1 + n1051ac2 + n1051ac3 + n1051ac4 + 
	n1051ac5 + n1051ac6 + 200.0*-1;
subject to acmiles:
	0 <= 2.67711*n1001ac1 + 2.67711*n1001ac2 + 2.67711*n1001ac3 + 2.80333*n1002ac1 
	+ 2.80333*n1002ac2 + 2.80333*n1002ac3 + 2.55338*n1003ac1 + 2.55338*n1003ac2 + 
	2.55338*n1003ac3 + 2.8919*n1004ac1 + 2.8919*n1004ac2 + 2.8919*n1004ac3 + 
	2.88537*n1005ac3 + 2.89192*n1105ac3 + 2.96093*n1006ac3 + 5.16595*n1007ac1 + 
	5.16595*n1007ac2 + 5.16595*n1007ac3 + 1.14422*n1008ac1 + 1.14422*n1008ac2 + 
	1.14422*n1008ac3 + 1.14422*n1008ac4 + 1.14422*n1008ac5 + 1.14422*n1008ac6 + 
	1.08196*n1009ac1 + 1.08196*n1009ac2 + 1.08196*n1009ac3 + 1.08196*n1009ac4 + 
	1.08196*n1009ac5 + 0.12622*n1010ac1 + 0.12622*n1010ac2 + 0.12622*n1010ac3 + 
	0.12622*n1010ac4 + 0.12622*n1010ac5 + 0.12622*n1010ac6 + 1.01801*n1011ac1 + 
	1.01801*n1011ac2 + 1.01801*n1011ac3 + 1.01801*n1011ac4 + 1.01801*n1011ac5 + 
	1.01801*n1011ac6 + 0.95575*n1012ac1 + 0.95575*n1012ac2 + 0.95575*n1012ac3 + 
	0.95575*n1012ac4 + 0.95575*n1012ac5 + 1.03338*n1013ac3 + 1.03338*n1013ac4 + 
	1.03338*n1013ac5 + 1.03338*n1013ac6 + 0.9975*n1014ac3 + 0.9975*n1014ac4 + 
	0.9975*n1014ac5 + 0.9975*n1014ac6 + 1.04277*n1015ac3 + 1.04277*n1015ac4 + 
	1.04277*n1015ac5 + 1.04277*n1015ac6 + 1.00605*n1016ac3 + 1.00605*n1016ac4 + 
	1.00605*n1016ac5 + 1.00605*n1016ac6 + 1.00997*n1017ac3 + 1.00997*n1017ac4 + 
	1.00997*n1017ac5 + 1.00997*n1017ac6 + 0.68996*n1018ac1 + 0.68996*n1018ac2 + 
	0.68996*n1018ac3 + 0.68996*n1018ac4 + 0.68996*n1018ac5 + 0.68996*n1018ac6 + 
	0.99393*n1019ac1 + 0.99393*n1019ac2 + 0.99393*n1019ac3 + 0.99393*n1019ac4 + 
	0.99393*n1019ac5 + 0.56374*n1020ac1 + 0.56374*n1020ac2 + 0.56374*n1020ac3 + 
	0.56374*n1020ac4 + 0.56374*n1020ac5 + 0.56374*n1020ac6 + 0.86771*n1021ac1 + 
	0.86771*n1021ac2 + 0.86771*n1021ac3 + 0.86771*n1021ac4 + 0.86771*n1021ac5 + 
	7.64184*n1022ac1 + 5.13904*n1023ac1 + 5.96857*n1026ac1 + 7.05137*n1027ac1 + 
	6.77512*n1028ac1 + 7.14547*n1029ac1 + 7.11364*n1030ac1 + 2.78679*n1032ac1 + 
	2.78679*n1032ac2 + 2.78679*n1032ac3 + 2.78679*n1032ac4 + 2.78679*n1032ac5 + 
	2.5415*n1033ac1 + 2.5415*n1033ac2 + 2.5415*n1033ac3 + 2.5415*n1033ac4 + 
	2.5415*n1033ac5 + 2.39258*n1034ac1 + 2.39258*n1034ac2 + 2.39258*n1034ac3 + 
	2.09179*n1035ac1 + 2.09179*n1035ac2 + 2.09179*n1035ac3 + 2.09179*n1035ac4 + 
	2.09179*n1035ac5 + 2.54144*n1036ac1 + 2.54144*n1036ac2 + 2.54144*n1036ac3 + 
	1.24837*n1037ac4 + 1.24837*n1037ac5 + 0.93376*n1038ac4 + 0.93376*n1038ac5 + 
	1.1293*n1039ac4 + 1.1293*n1039ac5 + 0.56879*n1040ac4 + 0.56879*n1040ac5 + 
	0.56879*n1040ac6 + 0.31461*n1041ac4 + 0.31461*n1041ac5 + 0.31461*n1041ac6 + 
	0.25418*n1042ac4 + 0.25418*n1042ac5 + 0.25418*n1042ac6 + 2.48884*n1043ac1 + 
	2.48884*n1043ac2 + 2.48884*n1043ac3 + 2.69845*n1044ac1 + 2.69845*n1044ac2 + 
	2.69845*n1044ac3 + 3.02542*n1046ac3 + 5.59035*n1047ac1 + 5.59035*n1047ac2 + 
	5.59035*n1047ac3 + 1.00139*n1050ac3 + 1.00139*n1050ac4 + 1.00139*n1050ac5 + 
	0.33852*n1051ac1 + 0.33852*n1051ac2 + 0.33852*n1051ac3 + 0.33852*n1051ac4 + 
	0.33852*n1051ac5 + 0.33852*n1051ac6;
subject to asmiles:
	0 <= 1086.90576*n1001ac1 + 690.69434*n1001ac2 + 492.58813*n1001ac3 + 
	1138.1499*n1002ac1 + 723.25806*n1002ac2 + 515.81177*n1002ac3 + 
	1036.67114*n1003ac1 + 658.77124*n1003ac2 + 469.82129*n1003ac3 + 
	1174.11182*n1004ac1 + 746.11035*n1004ac2 + 532.10986*n1004ac3 + 
	530.90723*n1005ac3 + 532.11377*n1105ac3 + 544.81104*n1006ac3 + 
	2097.37598*n1007ac1 + 1332.81494*n1007ac2 + 950.53516*n1007ac3 + 
	464.55518*n1008ac1 + 295.20972*n1008ac2 + 210.53731*n1008ac3 + 
	109.84555*n1008ac4 + 141.8838*n1008ac5 + 108.70131*n1008ac6 + 
	439.27612*n1009ac1 + 279.14575*n1009ac2 + 199.08081*n1009ac3 + 
	103.86824*n1009ac4 + 134.16312*n1009ac5 + 51.24359*n1010ac1 + 32.56366*n1010ac2 
	+ 23.22369*n1010ac3 + 12.11671*n1010ac4 + 15.65075*n1010ac5 + 11.9905*n1010ac6 
	+ 413.31152*n1011ac1 + 262.64624*n1011ac2 + 187.31361*n1011ac3 + 
	97.72881*n1011ac4 + 126.23305*n1011ac5 + 96.7108*n1011ac6 + 388.03247*n1012ac1 
	+ 246.58224*n1012ac2 + 175.85712*n1012ac3 + 91.7515*n1012ac4 + 
	118.51237*n1012ac5 + 190.14275*n1013ac3 + 99.20493*n1013ac4 + 
	128.13968*n1013ac5 + 98.17149*n1013ac6 + 183.53999*n1014ac3 + 95.75999*n1014ac4 
	+ 123.68999*n1014ac5 + 94.7625*n1014ac6 + 191.86919*n1015ac3 + 
	100.10562*n1015ac4 + 129.30312*n1015ac5 + 99.06287*n1015ac6 + 
	185.11406*n1016ac3 + 96.58124*n1016ac4 + 124.75075*n1016ac5 + 95.57518*n1016ac6 
	+ 185.8343*n1017ac3 + 96.95699*n1017ac4 + 125.23611*n1017ac5 + 
	95.94705*n1017ac6 + 280.12378*n1018ac1 + 178.00974*n1018ac2 + 
	126.95268*n1018ac3 + 66.23618*n1018ac4 + 85.55505*n1018ac5 + 65.54617*n1018ac6 
	+ 403.53516*n1019ac1 + 256.43359*n1019ac2 + 182.88293*n1019ac3 + 
	95.41718*n1019ac4 + 123.24718*n1019ac5 + 228.88037*n1020ac1 + 
	145.44612*n1020ac2 + 103.72899*n1020ac3 + 54.11951*n1020ac4 + 69.90431*n1020ac5 
	+ 53.55576*n1020ac6 + 352.2915*n1021ac1 + 223.87006*n1021ac2 + 
	159.65924*n1021ac3 + 83.30049*n1021ac4 + 107.59644*n1021ac5 + 
	3102.58496*n1022ac1 + 2086.44897*n1023ac1 + 2423.23682*n1026ac1 + 
	2862.85791*n1027ac1 + 2750.69678*n1028ac1 + 2901.05981*n1029ac1 + 
	2888.13794*n1030ac1 + 1131.43799*n1032ac1 + 718.99268*n1032ac2 + 
	512.76978*n1032ac3 + 267.53198*n1032ac4 + 345.56226*n1032ac5 + 
	1031.8501*n1033ac1 + 655.70752*n1033ac2 + 467.63647*n1033ac3 + 
	243.98424*n1033ac4 + 315.14624*n1033ac5 + 971.38599*n1034ac1 + 
	617.28467*n1034ac2 + 440.23389*n1034ac3 + 849.26563*n1035ac1 + 
	539.68115*n1035ac2 + 384.88867*n1035ac3 + 200.81155*n1035ac4 + 
	259.38159*n1035ac5 + 1031.82446*n1036ac1 + 655.69141*n1036ac2 + 
	467.62476*n1036ac3 + 119.84406*n1037ac4 + 154.79855*n1037ac5 + 
	89.64105*n1038ac4 + 115.78636*n1038ac5 + 108.41281*n1039ac4 + 
	140.03325*n1039ac5 + 54.60388*n1040ac4 + 70.53*n1040ac5 + 54.0351*n1040ac6 + 
	30.20299*n1041ac4 + 39.01219*n1041ac5 + 29.88837*n1041ac6 + 24.40089*n1042ac4 + 
	31.51784*n1042ac5 + 24.14673*n1042ac6 + 1010.47021*n1043ac1 + 
	642.12134*n1043ac2 + 457.94702*n1043ac3 + 1095.56982*n1044ac1 + 
	696.19995*n1044ac2 + 496.51465*n1044ac3 + 556.67676*n1046ac3 + 
	2269.68188*n1047ac1 + 1442.30981*n1047ac2 + 1028.62451*n1047ac3 + 
	184.25568*n1050ac3 + 96.13336*n1050ac4 + 124.1723*n1050ac5 + 137.44067*n1051ac1 
	+ 87.33911*n1051ac2 + 62.28839*n1051ac3 + 32.49829*n1051ac4 + 41.97696*n1051ac5 
	+ 32.15976*n1051ac6;
subject to passngrs:
	0 <= pboshnl0 + pboshnl1 - pboshnl7 - pboshnl8 + pboslax0 + pboslax1 - pboslax7 
	+ pbossea0 + pbossea1 + pbossea2 + pbossfo0 + pbossfo1 - pbostpe1 - pbostpe2 - 
	pbostyo1 - pbostyo2 + pbosyul0 + pbosyul1 + pbosyul2 + pbosyul3 + pbosyul4 + 
	pbosyvr0 + pbosyvr1 + pbosyvr2 + pbosywg0 + pbosywg1 - pbosywg7 + pbosyyz0 + 
	pbosyyz1 + pburoak0 + pburoak1 + pbursea0 + pbursea1 + pbursfo0 + phnllax0 + 
	phnllax1 + phnllax2 + phnllax3 + phnllon0 - phnllon6 + phnlpar0 - phnlpar6 + 
	phnlsea0 + phnlsea1 + phnlsea2 + phnlsfo0 + phnlsfo1 - phnlsfo7 + phnlyvr0 + 
	phnlyvr1 - phnlyvr7 + plassea0 + plassea1 + plasyvr0 - plasyvr6 + plaxoak0 + 
	plaxoak1 + plaxoak2 + plaxsea0 + plaxsea1 + plaxsea2 + plaxsea3 + plaxsea4 + 
	plaxsea5 + plaxsea6 + plaxsea7 + plaxsea8 + plaxsea9 + plaxsfo0 + plaxsfo1 + 
	plaxsfo2 + plaxsfo3 + plaxsfo4 + plaxsfo5 + plaxsfo6 + plaxtpe0 + plaxtpe1 + 
	plaxtpe2 - plaxtpe8 + plaxtyo0 + plaxtyo1 + plaxtyo2 - plaxtyo8 + plaxyvr0 + 
	plaxyvr1 + plaxyvr2 + plonpar0 + plonpar1 + plonpar2 + plonsea0 + plonsea1 + 
	plonyvr0 + plonyvr1 - plonyvr7 + poakont0 + poakont1 + poaksea0 + poaksea1 + 
	poaksea2 + pontsfo0 + pontsea0 + pontsea1 + pontsea2 + pparsea0 + pparsea1 + 
	pparyvr0 + pparyvr1 - pparyvr7 + prnosea0 + prnosea1 + prnoyvr0 + pseasfo0 + 
	pseasfo1 + pseasfo2 + pseasfo3 + pseasfo4 + pseasfo5 + pseasfo6 + pseatpe0 + 
	pseatpe1 + pseatpe2 + pseatpe3 + pseatyo0 + pseatyo1 + pseatyo2 + pseatyo3 + 
	pseayvr0 + pseayvr1 + pseayvr2 + pseayvr3 + pseayvr4 + pseayvr5 + pseayvr6 + 
	pseayvr7 + pseayvr8 + psfotpe0 + psfotpe1 + psfotpe2 - psfotpe8 + psfotyo0 + 
	psfotyo1 + psfotyo2 - psfotyo8 + psfoyvr0 + psfoyvr1 + ptpetyo0 + ptpetyo1 + 
	ptpetyo2 + ptpetyo3 + ptpeyvr0 + ptyoyvr0 + pyulyvr0 + pyulyvr1 + pyulyvr2 + 
	pyulyvr3 + pyulywg0 + pyulywg1 + pyulywg2 + pyulywg3 + pyulyyz0 + pyulyyz1 + 
	pyulyyz2 + pyulyyz3 + pyulyyz4 + pyvrywg0 + pyvrywg1 + pyvrywg2 + pyvryyz0 + 
	pyvryyz1 + pyvryyz2 + pywgyyz0 + pywgyyz1 + pywgyyz2 + pywgyyz3 + pbosoak0 - 
	pbosoak6 - pbosbur1 - 2.0*pbosbur2 - pbosont1 - 2.0*pbosont2 - pburyvr1 - 
	pburtyo1 - pburtpe1 + pburhnl0 - pburhnl6 + phnloak0 + phnloak1 + phnloak2 - 
	phnloak8 + phnlont0 - phnlont6 - phnlywg1 - phnlyyz1 - phnlyul1 - plastyo1 - 
	plastpe1 + plaxlon0 - plaxlon6 - plaxlon7 + plaxpar0 - plaxpar6 - plaxpar7 - 
	pburlon1 - pburpar1 - plonont1 - plonoak1 - poakpar1 - poaktyo1 - poaktpe1 - 
	pontpar1 - ponttyo1 - ponttpe1 - pparsfo1 - prnotyo1 - prnotpe1 - ptpeywg1 - 
	ptpeyyz1 - ptpeyul1 - ptyoyul1 - ptyoyyz1 - ptyoywg1 + plaxont0;
subject to rpmiles:
	0 <= 5.16595*pboshnl0 + 5.59035*pboshnl1 + 3.02542*pboslax0 + 3.03697*pboslax1 
	+ 2.48884*pbossea0 + 2.78679*pbossea1 + 2.48884*pbossea2 + 2.69845*pbossfo0 + 
	2.69845*pbossfo1 + 0.25418*pbosyul0 + 0.25418*pbosyul1 + 0.25418*pbosyul2 + 
	0.25418*pbosyul3 + 0.25418*pbosyul4 + 2.66058*pbosyvr0 + 2.5415*pbosyvr1 + 
	2.54144*pbosyvr2 + 1.50255*pbosywg0 + 1.38348*pbosywg1 + 0.56879*pbosyyz0 + 
	0.56879*pbosyyz1 + 0.32608*pburoak0 + 0.32608*pburoak1 + 0.9975*pbursea0 + 
	1.00605*pbursea1 + 0.32657*pbursfo0 + 2.55338*phnllax0 + 2.55338*phnllax1 + 
	2.55338*phnllax2 + 2.55338*phnllax3 + 7.41547*phnllon0 + 7.64184*phnlpar0 + 
	2.67711*phnlsea0 + 2.67711*phnlsea1 + 2.67711*phnlsea2 + 2.8919*phnlsfo0 + 
	2.8919*phnlsfo1 + 2.80333*phnlyvr0 + 2.70531*phnlyvr1 + 0.86771*plassea0 + 
	0.86771*plassea1 + 0.99393*plasyvr0 + 0.33855*plaxoak0 + 0.33855*plaxoak1 + 
	0.33855*plaxoak2 + 1.01801*plaxsea0 + 0.95575*plaxsea1 + 1.01801*plaxsea2 + 
	0.95575*plaxsea3 + 1.00997*plaxsea4 + 0.95575*plaxsea5 + 0.95575*plaxsea6 + 
	1.01801*plaxsea7 + 1.01801*plaxsea8 + 0.95575*plaxsea9 + 0.33852*plaxsfo0 + 
	0.33852*plaxsfo1 + 0.33852*plaxsfo2 + 0.33852*plaxsfo3 + 0.33852*plaxsfo4 + 
	0.33852*plaxsfo5 + 0.33852*plaxsfo6 + 7.05137*plaxtpe0 + 7.14547*plaxtpe1 + 
	7.11364*plaxtpe2 + 5.74891*plaxtyo0 + 5.843*plaxtyo1 + 5.81118*plaxtyo2 + 
	1.14422*plaxyvr0 + 1.08196*plaxyvr1 + 1.14422*plaxyvr2 + 0.22636*plonpar0 + 
	0.22636*plonpar1 + 0.22636*plonpar2 + 4.78646*plonsea0 + 4.78646*plonsea1 + 
	4.71017*plonyvr0 + 4.91267*plonyvr1 + 0.36196*poakont0 + 0.36196*poakont1 + 
	0.67142*poaksea0 + 0.67142*poaksea1 + 0.67142*poaksea2 + 0.36328*pontsfo0 + 
	1.03338*pontsea0 + 1.04277*pontsea1 + 1.00139*pontsea2 + 5.01282*pparsea0 + 
	5.01282*pparsea1 + 4.93653*pparyvr0 + 5.13904*pparyvr1 + 0.56374*prnosea0 + 
	0.56374*prnosea1 + 0.68996*prnoyvr0 + 0.67949*pseasfo0 + 0.67949*pseasfo1 + 
	0.67949*pseasfo2 + 0.67949*pseasfo3 + 0.67949*pseasfo4 + 0.67949*pseasfo5 + 
	0.67949*pseasfo6 + 6.09563*pseatpe0 + 6.09563*pseatpe1 + 6.12746*pseatpe2 + 
	6.09563*pseatpe3 + 4.79317*pseatyo0 + 4.79317*pseatyo1 + 4.825*pseatyo2 + 
	4.79317*pseatyo3 + 0.12622*pseayvr0 + 0.12622*pseayvr1 + 0.12622*pseayvr2 + 
	0.12622*pseayvr3 + 0.12622*pseayvr4 + 0.12622*pseayvr5 + 0.12622*pseayvr6 + 
	0.12622*pseayvr7 + 0.12622*pseayvr8 + 6.77512*psfotpe0 + 6.80694*psfotpe1 + 
	6.77512*psfotpe2 + 5.47265*psfotyo0 + 5.50448*psfotyo1 + 5.47265*psfotyo2 + 
	0.8057*psfoyvr0 + 0.8057*psfoyvr1 + 1.30247*ptpetyo0 + 1.30247*ptpetyo1 + 
	1.30247*ptpetyo2 + 1.30247*ptpetyo3 + 6.00125*ptpeyvr0 + 4.69878*ptyoyvr0 + 
	2.4064*pyulyvr0 + 2.28733*pyulyvr1 + 2.39258*pyulyvr2 + 2.28726*pyulyvr3 + 
	1.24837*pyulywg0 + 1.1293*pyulywg1 + 1.24837*pyulywg2 + 1.1293*pyulywg3 + 
	0.31461*pyulyyz0 + 0.31461*pyulyyz1 + 0.31461*pyulyyz2 + 0.31461*pyulyyz3 + 
	0.31461*pyulyyz4 + 1.15803*pyvrywg0 + 1.15803*pyvrywg1 + 1.15803*pyvrywg2 + 
	2.09179*pyvryyz0 + 2.07796*pyvryyz1 + 2.09179*pyvryyz2 + 0.93376*pywgyyz0 + 
	0.93376*pywgyyz1 + 0.93376*pywgyyz2 + 0.93376*pywgyyz3 + 2.68687*pbosoak0 + 
	2.55929*pburhnl0 + 2.88537*phnloak0 + 2.89192*phnloak1 + 2.96093*phnloak2 + 
	2.59897*phnlont0 + 5.7422*plaxlon0 + 5.96857*plaxpar0 + 0.04564*plaxont0;
subject to lfrpmasm:
	0 <= -5.16595*pboshnl0 - 5.59035*pboshnl1 - 3.02542*pboslax0 - 3.03697*pboslax1 
	- 2.48884*pbossea0 - 2.78679*pbossea1 - 2.48884*pbossea2 - 2.69845*pbossfo0 - 
	2.69845*pbossfo1 - 0.25418*pbosyul0 - 0.25418*pbosyul1 - 0.25418*pbosyul2 - 
	0.25418*pbosyul3 - 0.25418*pbosyul4 - 2.66058*pbosyvr0 - 2.5415*pbosyvr1 - 
	2.54144*pbosyvr2 - 1.50255*pbosywg0 - 1.38348*pbosywg1 - 0.56879*pbosyyz0 - 
	0.56879*pbosyyz1 - 0.32608*pburoak0 - 0.32608*pburoak1 - 0.9975*pbursea0 - 
	1.00605*pbursea1 - 0.32657*pbursfo0 - 2.55338*phnllax0 - 2.55338*phnllax1 - 
	2.55338*phnllax2 - 2.55338*phnllax3 - 7.41547*phnllon0 - 7.64184*phnlpar0 - 
	2.67711*phnlsea0 - 2.67711*phnlsea1 - 2.67711*phnlsea2 - 2.8919*phnlsfo0 - 
	2.8919*phnlsfo1 - 2.80333*phnlyvr0 - 2.70531*phnlyvr1 - 0.86771*plassea0 - 
	0.86771*plassea1 - 0.99393*plasyvr0 - 0.33855*plaxoak0 - 0.33855*plaxoak1 - 
	0.33855*plaxoak2 - 1.01801*plaxsea0 - 0.95575*plaxsea1 - 1.01801*plaxsea2 - 
	0.95575*plaxsea3 - 1.00997*plaxsea4 - 0.95575*plaxsea5 - 0.95575*plaxsea6 - 
	1.01801*plaxsea7 - 1.01801*plaxsea8 - 0.95575*plaxsea9 - 0.33852*plaxsfo0 - 
	0.33852*plaxsfo1 - 0.33852*plaxsfo2 - 0.33852*plaxsfo3 - 0.33852*plaxsfo4 - 
	0.33852*plaxsfo5 - 0.33852*plaxsfo6 - 7.05137*plaxtpe0 - 7.14547*plaxtpe1 - 
	7.11364*plaxtpe2 - 5.74891*plaxtyo0 - 5.843*plaxtyo1 - 5.81118*plaxtyo2 - 
	1.14422*plaxyvr0 - 1.08196*plaxyvr1 - 1.14422*plaxyvr2 - 0.22636*plonpar0 - 
	0.22636*plonpar1 - 0.22636*plonpar2 - 4.78646*plonsea0 - 4.78646*plonsea1 - 
	4.71017*plonyvr0 - 4.91267*plonyvr1 - 0.36196*poakont0 - 0.36196*poakont1 - 
	0.67142*poaksea0 - 0.67142*poaksea1 - 0.67142*poaksea2 - 0.36328*pontsfo0 - 
	1.03338*pontsea0 - 1.04277*pontsea1 - 1.00139*pontsea2 - 5.01282*pparsea0 - 
	5.01282*pparsea1 - 4.93653*pparyvr0 - 5.13904*pparyvr1 - 0.56374*prnosea0 - 
	0.56374*prnosea1 - 0.68996*prnoyvr0 - 0.67949*pseasfo0 - 0.67949*pseasfo1 - 
	0.67949*pseasfo2 - 0.67949*pseasfo3 - 0.67949*pseasfo4 - 0.67949*pseasfo5 - 
	0.67949*pseasfo6 - 6.09563*pseatpe0 - 6.09563*pseatpe1 - 6.12746*pseatpe2 - 
	6.09563*pseatpe3 - 4.79317*pseatyo0 - 4.79317*pseatyo1 - 4.825*pseatyo2 - 
	4.79317*pseatyo3 - 0.12622*pseayvr0 - 0.12622*pseayvr1 - 0.12622*pseayvr2 - 
	0.12622*pseayvr3 - 0.12622*pseayvr4 - 0.12622*pseayvr5 - 0.12622*pseayvr6 - 
	0.12622*pseayvr7 - 0.12622*pseayvr8 - 6.77512*psfotpe0 - 6.80694*psfotpe1 - 
	6.77512*psfotpe2 - 5.47265*psfotyo0 - 5.50448*psfotyo1 - 5.47265*psfotyo2 - 
	0.8057*psfoyvr0 - 0.8057*psfoyvr1 - 1.30247*ptpetyo0 - 1.30247*ptpetyo1 - 
	1.30247*ptpetyo2 - 1.30247*ptpetyo3 - 6.00125*ptpeyvr0 - 4.69878*ptyoyvr0 - 
	2.4064*pyulyvr0 - 2.28733*pyulyvr1 - 2.39258*pyulyvr2 - 2.28726*pyulyvr3 - 
	1.24837*pyulywg0 - 1.1293*pyulywg1 - 1.24837*pyulywg2 - 1.1293*pyulywg3 - 
	0.31461*pyulyyz0 - 0.31461*pyulyyz1 - 0.31461*pyulyyz2 - 0.31461*pyulyyz3 - 
	0.31461*pyulyyz4 - 1.15803*pyvrywg0 - 1.15803*pyvrywg1 - 1.15803*pyvrywg2 - 
	2.09179*pyvryyz0 - 2.07796*pyvryyz1 - 2.09179*pyvryyz2 - 0.93376*pywgyyz0 - 
	0.93376*pywgyyz1 - 0.93376*pywgyyz2 - 0.93376*pywgyyz3 - 2.68687*pbosoak0 - 
	2.55929*pburhnl0 - 2.88537*phnloak0 - 2.89192*phnloak1 - 2.96093*phnloak2 - 
	2.59897*phnlont0 - 5.7422*plaxlon0 - 5.96857*plaxpar0 - 0.04564*plaxont0 + 
	760.834032*n1001ac1 + 483.486038*n1001ac2 + 344.811691*n1001ac3 + 
	796.70493*n1002ac1 + 506.280642*n1002ac2 + 361.068239*n1002ac3 + 
	725.669798*n1003ac1 + 461.139868*n1003ac2 + 328.874903*n1003ac3 + 
	821.878274*n1004ac1 + 522.277245*n1004ac2 + 372.476902*n1004ac3 + 
	371.635061*n1005ac3 + 372.479639*n1105ac3 + 381.367728*n1006ac3 + 
	1468.163186*n1007ac1 + 932.970458*n1007ac2 + 665.374612*n1007ac3 + 
	325.188626*n1008ac1 + 206.646804*n1008ac2 + 147.376117*n1008ac3 + 
	76.891885*n1008ac4 + 99.31866*n1008ac5 + 76.090917*n1008ac6 + 
	307.493284*n1009ac1 + 195.402025*n1009ac2 + 139.356567*n1009ac3 + 
	72.707768*n1009ac4 + 93.914184*n1009ac5 + 35.870513*n1010ac1 + 
	22.794562*n1010ac2 + 16.256583*n1010ac3 + 8.481697*n1010ac4 + 
	10.955525*n1010ac5 + 8.39335*n1010ac6 + 289.318064*n1011ac1 + 
	183.852368*n1011ac2 + 131.119527*n1011ac3 + 68.410167*n1011ac4 + 
	88.363135*n1011ac5 + 67.69756*n1011ac6 + 271.622729*n1012ac1 + 
	172.607568*n1012ac2 + 123.099984*n1012ac3 + 64.22605*n1012ac4 + 
	82.958659*n1012ac5 + 133.099925*n1013ac3 + 69.443451*n1013ac4 + 
	89.697776*n1013ac5 + 68.720043*n1013ac6 + 128.477993*n1014ac3 + 
	67.031993*n1014ac4 + 86.582993*n1014ac5 + 66.33375*n1014ac6 + 
	134.308433*n1015ac3 + 70.073934*n1015ac4 + 90.512184*n1015ac5 + 
	69.344009*n1015ac6 + 129.579842*n1016ac3 + 67.606868*n1016ac4 + 
	87.325525*n1016ac5 + 66.902626*n1016ac6 + 130.08401*n1017ac3 + 
	67.869893*n1017ac4 + 87.665277*n1017ac5 + 67.162935*n1017ac6 + 
	196.086646*n1018ac1 + 124.606818*n1018ac2 + 88.866876*n1018ac3 + 
	46.365326*n1018ac4 + 59.888535*n1018ac5 + 45.882319*n1018ac6 + 
	282.474612*n1019ac1 + 179.503513*n1019ac2 + 128.018051*n1019ac3 + 
	66.792026*n1019ac4 + 86.273026*n1019ac5 + 160.216259*n1020ac1 + 
	101.812284*n1020ac2 + 72.610293*n1020ac3 + 37.883657*n1020ac4 + 
	48.933017*n1020ac5 + 37.489032*n1020ac6 + 246.60405*n1021ac1 + 
	156.709042*n1021ac2 + 111.761468*n1021ac3 + 58.310343*n1021ac4 + 
	75.317508*n1021ac5 + 2171.809472*n1022ac1 + 1460.514279*n1023ac1 + 
	1696.265774*n1026ac1 + 2004.000537*n1027ac1 + 1925.487746*n1028ac1 + 
	2030.741867*n1029ac1 + 2021.696558*n1030ac1 + 792.006593*n1032ac1 + 
	503.294876*n1032ac2 + 358.938846*n1032ac3 + 187.272386*n1032ac4 + 
	241.893582*n1032ac5 + 722.29507*n1033ac1 + 458.995264*n1033ac2 + 
	327.345529*n1033ac3 + 170.788968*n1033ac4 + 220.602368*n1033ac5 + 
	679.970193*n1034ac1 + 432.099269*n1034ac2 + 308.163723*n1034ac3 + 
	594.485941*n1035ac1 + 377.776805*n1035ac2 + 269.422069*n1035ac3 + 
	140.568085*n1035ac4 + 181.567113*n1035ac5 + 722.277122*n1036ac1 + 
	458.983987*n1036ac2 + 327.337332*n1036ac3 + 83.890842*n1037ac4 + 
	108.358985*n1037ac5 + 62.748735*n1038ac4 + 81.050452*n1038ac5 + 
	75.888967*n1039ac4 + 98.023275*n1039ac5 + 38.222716*n1040ac4 + 49.371*n1040ac5 
	+ 37.82457*n1040ac6 + 21.142093*n1041ac4 + 27.308533*n1041ac5 + 
	20.921859*n1041ac6 + 17.080623*n1042ac4 + 22.062488*n1042ac5 + 
	16.902711*n1042ac6 + 707.329147*n1043ac1 + 449.484938*n1043ac2 + 
	320.562914*n1043ac3 + 766.898874*n1044ac1 + 487.339965*n1044ac2 + 
	347.560255*n1044ac3 + 389.673732*n1046ac3 + 1588.777316*n1047ac1 + 
	1009.616867*n1047ac2 + 720.037157*n1047ac3 + 128.978976*n1050ac3 + 
	67.293352*n1050ac4 + 86.92061*n1050ac5 + 96.208469*n1051ac1 + 
	61.137377*n1051ac2 + 43.601873*n1051ac3 + 22.748803*n1051ac4 + 
	29.383872*n1051ac5 + 22.511832*n1051ac6;
subject to flav1:
	grdtimo1 + grdtimn1 + 5.53037*n1001ac1 + 6.41046*n1002ac1 + 5.30481*n1003ac1 + 
	6.57193*n1004ac1 + 10.71753*n1007ac1 + 4.03592*n1008ac1 + 3.27241*n1009ac1 + 
	0.88009*n1010ac1 + 3.15583*n1011ac1 + 2.39232*n1012ac1 + 2.5578*n1018ac1 + 
	3.11193*n1019ac1 + 1.67771*n1020ac1 + 2.23184*n1021ac1 + 15.88107*n1022ac1 + 
	11.31847*n1023ac1 + 12.8307*n1026ac1 + 14.80466*n1027ac1 + 14.30104*n1028ac1 + 
	16.27618*n1029ac1 + 15.56816*n1030ac1 + 8.33032*n1032ac1 + 6.58316*n1033ac1 + 
	5.66166*n1034ac1 + 5.11333*n1035ac1 + 5.93304*n1036ac1 + 5.18716*n1043ac1 + 
	5.56927*n1044ac1 + 12.14121*n1047ac1 + 1.26713*n1051ac1 + 10.5*-1 = 0;
subject to flav2:
	grdtimo2 + grdtimn2 + 5.5491*n1001ac2 + 6.35099*n1002ac2 + 5.31883*n1003ac2 + 
	6.51583*n1004ac2 + 10.74784*n1007ac2 + 3.8304*n1008ac2 + 3.14753*n1009ac2 + 
	0.80189*n1010ac2 + 3.02851*n1011ac2 + 2.34564*n1012ac2 + 2.41802*n1018ac2 + 
	2.9837*n1019ac2 + 1.61613*n1020ac2 + 2.18181*n1021ac2 + 8.02122*n1032ac2 + 
	6.43074*n1033ac2 + 5.58658*n1034ac2 + 5.02682*n1035ac2 + 5.86362*n1036ac2 + 
	5.19874*n1043ac2 + 5.58881*n1044ac2 + 12.10464*n1047ac2 + 1.19699*n1051ac2 + 
	13.65*-1 = 0;
subject to flav3:
	grdtimo3 + grdtimn3 + 5.74745*n1001ac3 + 6.41862*n1002ac3 + 5.50122*n1003ac3 + 
	6.59488*n1004ac3 + 6.58188*n1005ac3 + 6.59493*n1105ac3 + 6.73225*n1006ac3 + 
	11.12024*n1007ac3 + 3.53701*n1008ac3 + 2.9931*n1009ac3 + 0.67117*n1010ac3 + 
	2.86584*n1011ac3 + 2.32193*n1012ac3 + 2.89643*n1013ac3 + 2.82502*n1014ac3 + 
	2.91511*n1015ac3 + 2.84205*n1016ac3 + 2.84984*n1017ac3 + 2.21302*n1018ac3 + 
	2.81792*n1019ac3 + 1.54185*n1020ac3 + 2.14675*n1021ac3 + 7.64572*n1032ac3 + 
	6.31759*n1033ac3 + 5.60122*n1034ac3 + 5.00266*n1035ac3 + 5.89746*n1036ac3 + 
	5.3728*n1043ac3 + 5.78991*n1044ac3 + 6.86058*n1046ac3 + 12.3848*n1047ac3 + 
	2.83276*n1050ac3 + 1.09366*n1051ac3 + 23.5*-1 = 0;
subject to flav4:
	grdtimo4 + grdtimn4 + 3.40841*n1008ac4 + 2.90163*n1009ac4 + 0.6329*n1010ac4 + 
	2.77551*n1011ac4 + 2.26873*n1012ac4 + 2.80583*n1013ac4 + 2.73507*n1014ac4 + 
	2.82434*n1015ac4 + 2.75194*n1016ac4 + 2.75966*n1017ac4 + 2.1286*n1018ac4 + 
	2.72803*n1019ac4 + 1.4957*n1020ac4 + 2.09513*n1021ac4 + 7.41556*n1032ac4 + 
	6.16384*n1033ac4 + 4.893*n1035ac4 + 3.2298*n1037ac4 + 2.22538*n1038ac4 + 
	2.61098*n1039ac4 + 1.88965*n1040ac4 + 1.00442*n1041ac4 + 0.88524*n1042ac4 + 
	2.74274*n1050ac4 + 1.05157*n1051ac4 + 21.75*-1 = 0;
subject to flav5:
	grdtimo5 + grdtimn5 + 3.42334*n1008ac5 + 2.92913*n1009ac5 + 0.62384*n1010ac5 + 
	2.7995*n1011ac5 + 2.30529*n1012ac5 + 2.83067*n1013ac5 + 2.75793*n1014ac5 + 
	2.84969*n1015ac5 + 2.77527*n1016ac5 + 2.78321*n1017ac5 + 2.13455*n1018ac5 + 
	2.75069*n1019ac5 + 1.51071*n1020ac5 + 2.12685*n1021ac5 + 7.48883*n1032ac5 + 
	6.25562*n1033ac5 + 4.97605*n1035ac5 + 3.26646*n1037ac5 + 2.26073*n1038ac5 + 
	2.65709*n1039ac5 + 1.88894*n1040ac5 + 1.00572*n1041ac5 + 0.88321*n1042ac5 + 
	2.76582*n1050ac5 + 1.05419*n1051ac5 + 21.75*-1 = 0;
subject to flav6:
	grdtimo6 + grdtimn6 + 3.59515*n1008ac6 + 0.61067*n1010ac6 + 2.98447*n1011ac6 + 
	3.01988*n1013ac6 + 2.93724*n1014ac6 + 3.04149*n1015ac6 + 2.95694*n1016ac6 + 
	2.96596*n1017ac6 + 2.22898*n1018ac6 + 1.6183*n1020ac6 + 1.94992*n1040ac6 + 
	1.04456*n1041ac6 + 0.90537*n1042ac6 + 1.09962*n1051ac6 + 24.3*-1 = 0;
subject to lf1001s1:
	0 <= -phnlsea0 + 284.0*n1001ac1 + 180.0*n1001ac2 + 128.0*n1001ac3;
subject to lf1002s1:
	0 <= -phnlyvr0 - pseayvr0 + 243.0*n1002ac1 + 154.0*n1002ac2 + 110.0*n1002ac3;
subject to lf1002s2:
	0 <= -phnlsea1 - phnlyvr0 + 284.0*n1002ac1 + 180.0*n1002ac2 + 128.0*n1002ac3;
subject to lf1003s1:
	0 <= -phnllax0 + 243.0*n1003ac1 + 154.0*n1003ac2 + 110.0*n1003ac3;
subject to lf1004s1:
	0 <= -phnlsfo0 - plaxsfo0 + 243.0*n1004ac1 + 154.0*n1004ac2 + 110.0*n1004ac3;
subject to lf1004s2:
	0 <= -phnllax1 - phnlsfo0 + 243.0*n1004ac1 + 154.0*n1004ac2 + 110.0*n1004ac3;
subject to lf1005s1:
	0 <= -pburoak0 - phnloak0 + 128.0*n1005ac3;
subject to lf1005s2:
	0 <= -pburhnl0 - phnloak0 + 128.0*n1005ac3;
subject to lf1105s1:
	0 <= -plaxoak0 - phnloak1 + 128.0*n1105ac3;
subject to lf1105s2:
	0 <= -phnllax2 - phnloak1 + 110.0*n1105ac3;
subject to lf1006s1:
	0 <= -poakont0 - phnloak2 + 128.0*n1006ac3;
subject to lf1006s2:
	0 <= -phnloak2 - phnlont0 + 128.0*n1006ac3;
subject to lf1007s1:
	0 <= -pboshnl0 - pbossea0 + 284.0*n1007ac1 + 180.0*n1007ac2 + 128.0*n1007ac3;
subject to lf1007s2:
	0 <= -pboshnl0 - phnlsea2 + 284.0*n1007ac1 + 180.0*n1007ac2 + 128.0*n1007ac3;
subject to lf1008s1:
	0 <= -plaxyvr0 - pseayvr1 - psfoyvr0 + 243.0*n1008ac1 + 154.0*n1008ac2 + 
	110.0*n1008ac3 + 57.0*n1008ac4 + 74.0*n1008ac5 + 56.0*n1008ac6;
subject to lf1008s2:
	0 <= -plaxsea0 - plaxyvr0 - pseasfo0 - psfoyvr0 + 243.0*n1008ac1 + 
	154.0*n1008ac2 + 110.0*n1008ac3 + 57.0*n1008ac4 + 74.0*n1008ac5 + 56.0*n1008ac6;
subject to lf1008s3:
	0 <= -plaxsea0 - plaxsfo1 - plaxyvr0 + 243.0*n1008ac1 + 154.0*n1008ac2 + 
	110.0*n1008ac3 + 57.0*n1008ac4 + 74.0*n1008ac5 + 56.0*n1008ac6;
subject to lf1009s1:
	0 <= -plaxyvr1 - pseayvr2 + 243.0*n1009ac1 + 154.0*n1009ac2 + 110.0*n1009ac3 + 
	57.0*n1009ac4 + 74.0*n1009ac5;
subject to lf1009s2:
	0 <= -plaxsea1 - plaxyvr1 + 243.0*n1009ac1 + 154.0*n1009ac2 + 110.0*n1009ac3 + 
	57.0*n1009ac4 + 74.0*n1009ac5;
subject to lf1010s1:
	0 <= -pseayvr3 + 243.0*n1010ac1 + 154.0*n1010ac2 + 110.0*n1010ac3 + 
	57.0*n1010ac4 + 74.0*n1010ac5 + 56.0*n1010ac6;
subject to lf1011s1:
	0 <= -plaxsea2 - pseasfo1 + 243.0*n1011ac1 + 154.0*n1011ac2 + 110.0*n1011ac3 + 
	57.0*n1011ac4 + 74.0*n1011ac5 + 56.0*n1011ac6;
subject to lf1011s2:
	0 <= -plaxsea2 - plaxsfo2 + 243.0*n1011ac1 + 154.0*n1011ac2 + 110.0*n1011ac3 + 
	57.0*n1011ac4 + 74.0*n1011ac5 + 56.0*n1011ac6;
subject to lf1012s1:
	0 <= -plaxsea3 + 243.0*n1012ac1 + 154.0*n1012ac2 + 110.0*n1012ac3 + 
	57.0*n1012ac4 + 74.0*n1012ac5;
subject to lf1013s1:
	0 <= -poaksea0 - pontsea0 + 128.0*n1013ac3 + 67.0*n1013ac4 + 86.0*n1013ac5 + 
	66.0*n1013ac6;
subject to lf1013s2:
	0 <= -poakont1 - pontsea0 + 128.0*n1013ac3 + 67.0*n1013ac4 + 86.0*n1013ac5 + 
	66.0*n1013ac6;
subject to lf1014s1:
	0 <= -pbursea0 - poaksea1 + 128.0*n1014ac3 + 67.0*n1014ac4 + 86.0*n1014ac5 + 
	66.0*n1014ac6;
subject to lf1014s2:
	0 <= -pburoak1 - pbursea0 + 128.0*n1014ac3 + 67.0*n1014ac4 + 86.0*n1014ac5 + 
	66.0*n1014ac6;
subject to lf1015s1:
	0 <= -pontsea1 - pseasfo2 + 110.0*n1015ac3 + 57.0*n1015ac4 + 74.0*n1015ac5 + 
	56.0*n1015ac6;
subject to lf1015s2:
	0 <= -pontsfo0 - pontsea1 + 128.0*n1015ac3 + 67.0*n1015ac4 + 86.0*n1015ac5 + 
	66.0*n1015ac6;
subject to lf1016s1:
	0 <= -pbursea1 - pseasfo3 + 110.0*n1016ac3 + 57.0*n1016ac4 + 74.0*n1016ac5 + 
	56.0*n1016ac6;
subject to lf1016s2:
	0 <= -pbursea1 - pbursfo0 + 128.0*n1016ac3 + 67.0*n1016ac4 + 86.0*n1016ac5 + 
	66.0*n1016ac6;
subject to lf1017s1:
	0 <= -plaxsea4 - poaksea2 + 128.0*n1017ac3 + 67.0*n1017ac4 + 86.0*n1017ac5 + 
	66.0*n1017ac6;
subject to lf1017s2:
	0 <= -plaxoak1 - plaxsea4 + 128.0*n1017ac3 + 67.0*n1017ac4 + 86.0*n1017ac5 + 
	66.0*n1017ac6;
subject to lf1018s1:
	0 <= -prnoyvr0 - pseayvr4 + 243.0*n1018ac1 + 154.0*n1018ac2 + 110.0*n1018ac3 + 
	57.0*n1018ac4 + 74.0*n1018ac5 + 56.0*n1018ac6;
subject to lf1018s2:
	0 <= -prnosea0 - prnoyvr0 + 243.0*n1018ac1 + 154.0*n1018ac2 + 110.0*n1018ac3 + 
	57.0*n1018ac4 + 74.0*n1018ac5 + 56.0*n1018ac6;
subject to lf1019s1:
	0 <= -plasyvr0 - pseayvr5 + 243.0*n1019ac1 + 154.0*n1019ac2 + 110.0*n1019ac3 + 
	57.0*n1019ac4 + 74.0*n1019ac5;
subject to lf1019s2:
	0 <= -plassea0 - plasyvr0 + 243.0*n1019ac1 + 154.0*n1019ac2 + 110.0*n1019ac3 + 
	57.0*n1019ac4 + 74.0*n1019ac5;
subject to lf1020s1:
	0 <= -prnosea1 + 243.0*n1020ac1 + 154.0*n1020ac2 + 110.0*n1020ac3 + 
	57.0*n1020ac4 + 74.0*n1020ac5 + 56.0*n1020ac6;
subject to lf1021s1:
	0 <= -plassea1 + 243.0*n1021ac1 + 154.0*n1021ac2 + 110.0*n1021ac3 + 
	57.0*n1021ac4 + 74.0*n1021ac5;
subject to lf1022s1:
	0 <= -phnllon0 - phnlpar0 - phnlyvr1 + 284.0*n1022ac1;
subject to lf1022s2:
	0 <= -phnllon0 - phnlpar0 - plonyvr0 - pparyvr0 + 284.0*n1022ac1;
subject to lf1022s3:
	0 <= -phnlpar0 - plonpar0 - pparyvr0 + 284.0*n1022ac1;
subject to lf1023s1:
	0 <= -plonyvr1 - pparyvr1 - pseayvr6 + 243.0*n1023ac1;
subject to lf1023s2:
	0 <= -plonsea0 - plonyvr1 - pparsea0 - pparyvr1 + 243.0*n1023ac1;
subject to lf1023s3:
	0 <= -plonpar1 - pparsea0 - pparyvr1 + 284.0*n1023ac1;
subject to lf1026s1:
	0 <= -plaxsea5 - plaxlon0 - plaxpar0 + 243.0*n1026ac1;
subject to lf1026s2:
	0 <= -plonsea1 - pparsea1 - plaxlon0 - plaxpar0 + 243.0*n1026ac1;
subject to lf1026s3:
	0 <= -plonpar2 - pparsea1 - plaxpar0 + 284.0*n1026ac1;
subject to lf1027s1:
	0 <= -plaxsea6 - plaxtpe0 - plaxtyo0 + 243.0*n1027ac1;
subject to lf1027s2:
	0 <= -plaxtpe0 - plaxtyo0 - pseatpe0 - pseatyo0 + 263.0*n1027ac1;
subject to lf1027s3:
	0 <= -plaxtpe0 - pseatpe0 - ptpetyo0 + 284.0*n1027ac1;
subject to lf1028s1:
	0 <= -pseasfo4 - psfotpe0 - psfotyo0 + 243.0*n1028ac1;
subject to lf1028s2:
	0 <= -pseatpe1 - pseatyo1 - psfotpe0 - psfotyo0 + 263.0*n1028ac1;
subject to lf1028s3:
	0 <= -pseatpe1 - psfotpe0 - ptpetyo1 + 284.0*n1028ac1;
subject to lf1029s1:
	0 <= -plaxsea7 - plaxsfo3 - plaxtpe1 - plaxtyo1 - plaxyvr2 + 243.0*n1029ac1;
subject to lf1029s2:
	0 <= -plaxsea7 - plaxtpe1 - plaxtyo1 - plaxyvr2 - pseasfo5 - psfotpe1 - 
	psfotyo1 - psfoyvr1 + 243.0*n1029ac1;
subject to lf1029s3:
	0 <= -plaxtpe1 - plaxtyo1 - plaxyvr2 - pseatpe2 - pseatyo2 - pseayvr7 - 
	psfotpe1 - psfotyo1 - psfoyvr1 + 243.0*n1029ac1;
subject to lf1029s4:
	0 <= -plaxtpe1 - plaxtyo1 - pseatpe2 - pseatyo2 - psfotpe1 - psfotyo1 - 
	ptpeyvr0 - ptyoyvr0 + 304.0*n1029ac1;
subject to lf1029s5:
	0 <= -plaxtpe1 - pseatpe2 - psfotpe1 - ptpetyo2 - ptpeyvr0 + 284.0*n1029ac1;
subject to lf1030s1:
	0 <= -plaxsea8 - plaxsfo4 - plaxtpe2 - plaxtyo2 + 243.0*n1030ac1;
subject to lf1030s2:
	0 <= -plaxsea8 - plaxtpe2 - plaxtyo2 - pseasfo6 - psfotpe2 - psfotyo2 + 
	243.0*n1030ac1;
subject to lf1030s3:
	0 <= -plaxtpe2 - plaxtyo2 - pseatpe3 - pseatyo3 - psfotpe2 - psfotyo2 + 
	263.0*n1030ac1;
subject to lf1030s4:
	0 <= -plaxtpe2 - pseatpe3 - psfotpe2 - ptpetyo3 + 284.0*n1030ac1;
subject to lf1032s1:
	0 <= -pbossea1 - pseayvr8 + 243.0*n1032ac1 + 154.0*n1032ac2 + 110.0*n1032ac3 + 
	57.0*n1032ac4 + 74.0*n1032ac5;
subject to lf1032s2:
	0 <= -pbossea1 - pbosyvr0 - pyulyvr0 - pyvrywg0 - pyvryyz0 + 284.0*n1032ac1 + 
	180.0*n1032ac2 + 128.0*n1032ac3 + 67.0*n1032ac4 + 86.0*n1032ac5;
subject to lf1032s3:
	0 <= -pbossea1 - pbosyvr0 - pbosywg0 - pyulyvr0 - pyulywg0 - pyvryyz0 - 
	pywgyyz0 + 284.0*n1032ac1 + 180.0*n1032ac2 + 128.0*n1032ac3 + 67.0*n1032ac4 + 
	86.0*n1032ac5;
subject to lf1032s4:
	0 <= -pbossea1 - pbosyvr0 - pbosywg0 - pbosyyz0 - pyulyvr0 - pyulywg0 - 
	pyulyyz0 + 284.0*n1032ac1 + 180.0*n1032ac2 + 128.0*n1032ac3 + 67.0*n1032ac4 + 
	86.0*n1032ac5;
subject to lf1032s5:
	0 <= -pbossea1 - pbosyul0 - pbosyvr0 - pbosywg0 - pbosyyz0 + 243.0*n1032ac1 + 
	154.0*n1032ac2 + 110.0*n1032ac3 + 57.0*n1032ac4 + 74.0*n1032ac5;
subject to lf1033s1:
	0 <= -pbosyvr1 - pyulyvr1 - pyvrywg1 + 284.0*n1033ac1 + 180.0*n1033ac2 + 
	128.0*n1033ac3 + 67.0*n1033ac4 + 86.0*n1033ac5;
subject to lf1033s2:
	0 <= -pbosyvr1 - pbosywg1 - pyulyvr1 - pyulywg1 + 284.0*n1033ac1 + 
	180.0*n1033ac2 + 128.0*n1033ac3 + 67.0*n1033ac4 + 86.0*n1033ac5;
subject to lf1033s3:
	0 <= -pbosyul1 - pbosyvr1 - pbosywg1 + 243.0*n1033ac1 + 154.0*n1033ac2 + 
	110.0*n1033ac3 + 57.0*n1033ac4 + 74.0*n1033ac5;
subject to lf1034s1:
	0 <= -pyulyvr2 - pyvryyz1 + 284.0*n1034ac1 + 180.0*n1034ac2 + 128.0*n1034ac3;
subject to lf1034s2:
	0 <= -pyulyvr2 - pyulyyz1 + 284.0*n1034ac1 + 180.0*n1034ac2 + 128.0*n1034ac3;
subject to lf1035s1:
	0 <= -pyvrywg2 - pyvryyz2 + 284.0*n1035ac1 + 180.0*n1035ac2 + 128.0*n1035ac3 + 
	67.0*n1035ac4 + 86.0*n1035ac5;
subject to lf1035s2:
	0 <= -pyvryyz2 - pywgyyz1 + 284.0*n1035ac1 + 180.0*n1035ac2 + 128.0*n1035ac3 + 
	67.0*n1035ac4 + 86.0*n1035ac5;
subject to lf1036s1:
	0 <= -pbosyvr2 - pyulyvr3 + 284.0*n1036ac1 + 180.0*n1036ac2 + 128.0*n1036ac3;
subject to lf1036s2:
	0 <= -pbosyul2 - pbosyvr2 + 243.0*n1036ac1 + 154.0*n1036ac2 + 110.0*n1036ac3;
subject to lf1037s1:
	0 <= -pyulywg2 - pywgyyz2 + 67.0*n1037ac4 + 86.0*n1037ac5;
subject to lf1037s2:
	0 <= -pyulywg2 - pyulyyz2 + 67.0*n1037ac4 + 86.0*n1037ac5;
subject to lf1038s1:
	0 <= -pywgyyz3 + 67.0*n1038ac4 + 86.0*n1038ac5;
subject to lf1039s1:
	0 <= -pyulywg3 + 67.0*n1039ac4 + 86.0*n1039ac5;
subject to lf1040s1:
	0 <= -pbosyyz1 - pyulyyz3 + 67.0*n1040ac4 + 86.0*n1040ac5 + 66.0*n1040ac6;
subject to lf1040s2:
	0 <= -pbosyul3 - pbosyyz1 + 57.0*n1040ac4 + 74.0*n1040ac5 + 56.0*n1040ac6;
subject to lf1041s1:
	0 <= -pyulyyz4 + 67.0*n1041ac4 + 86.0*n1041ac5 + 66.0*n1041ac6;
subject to lf1042s1:
	0 <= -pbosyul4 + 57.0*n1042ac4 + 74.0*n1042ac5 + 56.0*n1042ac6;
subject to lf1043s1:
	0 <= -pbossea2 + 284.0*n1043ac1 + 180.0*n1043ac2 + 128.0*n1043ac3;
subject to lf1044s1:
	0 <= -pbossfo0 + 243.0*n1044ac1 + 154.0*n1044ac2 + 110.0*n1044ac3;
subject to lf1046s1:
	0 <= -pboslax0 - plaxoak2 + 128.0*n1046ac3;
subject to lf1046s2:
	0 <= -pboslax0 - pbosoak0 + 128.0*n1046ac3;
subject to lf1047s1:
	0 <= -pboshnl1 - phnllax3 - phnlsfo1 + 243.0*n1047ac1 + 154.0*n1047ac2 + 
	110.0*n1047ac3;
subject to lf1047s2:
	0 <= -pboshnl1 - pboslax1 - phnlsfo1 - plaxsfo5 + 243.0*n1047ac1 + 
	154.0*n1047ac2 + 110.0*n1047ac3;
subject to lf1047s3:
	0 <= -pboshnl1 - pboslax1 - pbossfo1 + 243.0*n1047ac1 + 154.0*n1047ac2 + 
	110.0*n1047ac3;
subject to lf1050s1:
	0 <= -plaxsea9 - pontsea2 + 110.0*n1050ac3 + 57.0*n1050ac4 + 74.0*n1050ac5;
subject to lf1050s2:
	0 <= -pontsea2 - plaxont0 + 128.0*n1050ac3 + 67.0*n1050ac4 + 86.0*n1050ac5;
subject to lf1051s1:
	0 <= -plaxsfo6 + 243.0*n1051ac1 + 154.0*n1051ac2 + 110.0*n1051ac3 + 
	57.0*n1051ac4 + 74.0*n1051ac5 + 56.0*n1051ac6;
subject to noptlon0:
	0 >= 2.0*n1022ac1 + 2.0*n1023ac1 + 2.0*n1026ac1 + 4.0*-1;
subject to noptlon1:
	0 <= 2.0*n1022ac1 + 2.0*n1023ac1 + 2.0*n1026ac1 + 2.0*-1;
subject to nopttyo0:
	0 >= 2.0*n1027ac1 + 2.0*n1028ac1 + 2.0*n1029ac1 + 2.0*n1030ac1 + 4.0*-1;
subject to nopttyo1:
	0 <= 2.0*n1027ac1 + 2.0*n1028ac1 + 2.0*n1029ac1 + 2.0*n1030ac1 + 2.0*-1;
subject to dmboshnl:
	0 <= pboshnl0 + pboshnl1 + pboshnl7 + pboshnl8 + 12.0*-1 <= 2.0;
subject to dmboslax:
	0 <= pboslax0 + pboslax1 + pboslax7 + 14.0*-1 <= 2.0;
subject to dmbossea:
	0 <= -pboshnl7 + pbossea0 + pbossea1 + pbossea2 - pbostpe1 - pbostyo1 + 
	45.0*-1 <= 5.0;
subject to dmbossfo:
	0 <= -pboslax7 + pbossfo0 + pbossfo1 - pbosoak6 - pbosbur1 - pbosbur2 - 
	pbosont1 - pbosont2 + 122.0*-1 <= 13.0;
subject to dmbostpe:
	0 <= pbostpe1 + pbostpe2 + -1 <= 1.0;
subject to dmbostyo:
	0 <= pbostyo1 + pbostyo2 + 3.0*-1 <= 3.0;
subject to dmbosyul:
	0 <= pbosyul0 + pbosyul1 + pbosyul2 + pbosyul3 + pbosyul4 - pbosywg7 + 
	676.0*-1 <= 68.0;
subject to dmbosyvr:
	0 <= -pboshnl8 - pbostpe2 - pbostyo2 + pbosyvr0 + pbosyvr1 + pbosyvr2 + 
	26.0*-1 <= 3.0;
subject to dmbosywg:
	0 <= pbosywg0 + pbosywg1 + pbosywg7 + 37.0*-1 <= 4.0;
subject to dmbosyyz:
	0 <= pbosyyz0 + pbosyyz1 + 215.0*-1 <= 22.0;
subject to dmburoak:
	0 <= pburoak0 + pburoak1 + 27.0*-1 <= 3.0;
subject to dmbursea:
	0 <= pbursea0 + pbursea1 - pburyvr1 - pburtyo1 - pburtpe1 - pburlon1 - pburpar1 
	+ 52.0*-1 <= 6.0;
subject to dmbursfo:
	0 <= pbursfo0 - pbosbur1 + 271.0*-1 <= 28.0;
subject to dmhnllax:
	0 <= phnllax0 + phnllax1 + phnllax2 + phnllax3 - phnlsfo7 - pburhnl6 - phnloak8 
	- phnlont6 + 297.0*-1 <= 30.0;
subject to dmhnllon:
	0 <= phnllon0 + phnllon6 + 5.0*-1 <= 5.0;
subject to dmhnlpar:
	0 <= phnlpar0 + phnlpar6 + -1 <= 1.0;
subject to dmhnlsea:
	0 <= -pboshnl7 - phnllon6 - phnlpar6 + phnlsea0 + phnlsea1 + phnlsea2 - 
	phnlyvr7 + 112.0*-1 <= 12.0;
subject to dmhnlsfo:
	0 <= phnlsfo0 + phnlsfo1 + phnlsfo7 + 35.0*-1 <= 4.0;
subject to dmhnlyvr:
	0 <= -pboshnl8 + phnlyvr0 + phnlyvr1 + phnlyvr7 - phnlywg1 - phnlyyz1 - 
	phnlyul1 + 67.0*-1 <= 7.0;
subject to dmlassea:
	0 <= plassea0 + plassea1 - plasyvr6 - plastyo1 - plastpe1 + 370.0*-1 <= 38.0;
subject to dmlasyvr:
	0 <= plasyvr0 + plasyvr6 + 37.0*-1 <= 4.0;
subject to dmlaxoak:
	0 <= plaxoak0 + plaxoak1 + plaxoak2 - phnloak8 + 78.0*-1 <= 8.0;
subject to dmlaxsea:
	0 <= plaxsea0 + plaxsea1 + plaxsea2 + plaxsea3 + plaxsea4 + plaxsea5 + plaxsea6 
	+ plaxsea7 + plaxsea8 + plaxsea9 - plaxtpe8 - plaxtyo8 - plaxlon6 - plaxpar6 + 
	813.0*-1 <= 82.0;
subject to dmlaxsfo:
	0 <= -pboslax7 - phnlsfo7 + plaxsfo0 + plaxsfo1 + plaxsfo2 + plaxsfo3 + 
	plaxsfo4 + plaxsfo5 + plaxsfo6 - pbosbur2 - pbosont2 + 2952.0*-1 <= 296.0;
subject to dmlaxtpe:
	0 <= plaxtpe0 + plaxtpe1 + plaxtpe2 + plaxtpe8 + 31.0*-1 <= 4.0;
subject to dmlaxtyo:
	0 <= plaxtyo0 + plaxtyo1 + plaxtyo2 + plaxtyo8 + 41.0*-1 <= 5.0;
subject to dmlaxyvr:
	0 <= plaxyvr0 + plaxyvr1 + plaxyvr2 - plaxlon7 - plaxpar7 + 193.0*-1 <= 20.0;
subject to dmlonpar:
	0 <= plonpar0 + plonpar1 + plonpar2 + 2.0*-1 <= 2.0;
subject to dmlonsea:
	0 <= -phnllon6 + plonsea0 + plonsea1 - plonyvr7 - plaxlon6 - pburlon1 - 
	plonont1 - plonoak1 + 92.0*-1 <= 10.0;
subject to dmlonyvr:
	0 <= plonyvr0 + plonyvr1 + plonyvr7 - plaxlon7 + 51.0*-1 <= 6.0;
subject to dmoakont:
	0 <= poakont0 + poakont1 + 13.0*-1 <= 2.0;
subject to dmoaksea:
	0 <= poaksea0 + poaksea1 + poaksea2 - plonoak1 - poakpar1 - poaktyo1 - poaktpe1 
	+ 110.0*-1 <= 12.0;
subject to dmontsfo:
	0 <= pontsfo0 - pbosont1 + 173.0*-1 <= 18.0;
subject to dmontsea:
	0 <= pontsea0 + pontsea1 + pontsea2 - plonont1 - pontpar1 - ponttyo1 - ponttpe1 
	+ 42.0*-1 <= 5.0;
subject to dmparsea:
	0 <= -phnlpar6 + pparsea0 + pparsea1 - pparyvr7 - plaxpar6 - pburpar1 - 
	poakpar1 - pontpar1 - pparsfo1 + 36.0*-1 <= 4.0;
subject to dmparyvr:
	0 <= pparyvr0 + pparyvr1 + pparyvr7 - plaxpar7 + 24.0*-1 <= 3.0;
subject to dmrnosea:
	0 <= prnosea0 + prnosea1 - prnotyo1 - prnotpe1 + 284.0*-1 <= 29.0;
subject to dmrnoyvr:
	0 <= prnoyvr0 + 67.0*-1 <= 7.0;
subject to dmseasfo:
	0 <= pseasfo0 + pseasfo1 + pseasfo2 + pseasfo3 + pseasfo4 + pseasfo5 + pseasfo6 
	- psfotpe8 - psfotyo8 - pparsfo1 + 1417.0*-1 <= 142.0;
subject to dmseatpe:
	0 <= -pbostpe1 - plaxtpe8 + pseatpe0 + pseatpe1 + pseatpe2 + pseatpe3 - 
	psfotpe8 - pburtpe1 - plastpe1 - poaktpe1 - ponttpe1 - prnotpe1 + 47.0*-1 <= 
	5.0;
subject to dmseatyo:
	0 <= -pbostyo1 - plaxtyo8 + pseatyo0 + pseatyo1 + pseatyo2 + pseatyo3 - 
	psfotyo8 - pburtyo1 - plastyo1 - poaktyo1 - ponttyo1 - prnotyo1 + 114.0*-1 <= 
	12.0;
subject to dmseayvr:
	0 <= -phnlyvr7 - plasyvr6 - plonyvr7 - pparyvr7 + pseayvr0 + pseayvr1 + 
	pseayvr2 + pseayvr3 + pseayvr4 + pseayvr5 + pseayvr6 + pseayvr7 + pseayvr8 - 
	pburyvr1 + 547.0*-1 <= 55.0;
subject to dmsfotpe:
	0 <= psfotpe0 + psfotpe1 + psfotpe2 + psfotpe8 + 6.0*-1 <= 6.0;
subject to dmsfotyo:
	0 <= psfotyo0 + psfotyo1 + psfotyo2 + psfotyo8 + 17.0*-1 <= 2.0;
subject to dmsfoyvr:
	0 <= psfoyvr0 + psfoyvr1 + 298.0*-1 <= 30.0;
subject to dmtpetyo:
	0 <= ptpetyo0 + ptpetyo1 + ptpetyo2 + ptpetyo3 + 111.0*-1 <= 12.0;
subject to dmtpeyvr:
	0 <= -pbostpe2 + ptpeyvr0 - ptpeywg1 - ptpeyyz1 - ptpeyul1 + 17.0*-1 <= 2.0;
subject to dmtyoyvr:
	0 <= -pbostyo2 + ptyoyvr0 - ptyoyul1 - ptyoyyz1 - ptyoywg1 + 42.0*-1 <= 5.0;
subject to dmyulyvr:
	0 <= pyulyvr0 + pyulyvr1 + pyulyvr2 + pyulyvr3 - phnlyul1 - ptpeyul1 - ptyoyul1 
	+ 262.0*-1 <= 27.0;
subject to dmyulywg:
	0 <= -pbosywg7 + pyulywg0 + pyulywg1 + pyulywg2 + pyulywg3 + 413.0*-1 <= 42.0;
subject to dmyulyyz:
	0 <= pyulyyz0 + pyulyyz1 + pyulyyz2 + pyulyyz3 + pyulyyz4 + 2612.0*-1 <= 
	262.0;
subject to dmyvrywg:
	0 <= pyvrywg0 + pyvrywg1 + pyvrywg2 - phnlywg1 - ptpeywg1 - ptyoywg1 + 
	375.0*-1 <= 38.0;
subject to dmyvryyz:
	0 <= pyvryyz0 + pyvryyz1 + pyvryyz2 - phnlyyz1 - ptpeyyz1 - ptyoyyz1 + 
	318.0*-1 <= 32.0;
subject to dmywgyyz:
	0 <= pywgyyz0 + pywgyyz1 + pywgyyz2 + pywgyyz3 + 278.0*-1 <= 28.0;
subject to dmbosoak:
	0 <= pbosoak0 + pbosoak6 + 11.0*-1 <= 2.0;
subject to dmbosbur:
	0 <= pbosbur1 + pbosbur2 + 7.0*-1 <= 7.0;
subject to dmbosont:
	0 <= pbosont1 + pbosont2 + 4.0*-1 <= 4.0;
subject to dmburyvr:
	0 <= pburyvr1 + 26.0*-1 <= 3.0;
subject to dmburtyo:
	0 <= pburtyo1 + 2.0*-1 <= 2.0;
subject to dmburtpe:
	0 <= pburtpe1 + -1 <= 1.0;
subject to dmburhnl:
	0 <= pburhnl0 + pburhnl6 + 11.0*-1 <= 2.0;
subject to dmhnloak:
	0 <= phnloak0 + phnloak1 + phnloak2 + phnloak8 + 24.0*-1 <= 3.0;
subject to dmhnlont:
	0 <= phnlont0 + phnlont6 + 16.0*-1 <= 2.0;
subject to dmhnlywg:
	0 <= phnlywg1 + 3.0*-1 <= 3.0;
subject to dmhnlyyz:
	0 <= phnlyyz1 + 24.0*-1 <= 3.0;
subject to dmhnlyul:
	0 <= phnlyul1 + 40.0*-1 <= 5.0;
subject to dmlastyo:
	0 <= plastyo1 + 5.0*-1 <= 5.0;
subject to dmlastpe:
	0 <= plastpe1 + -1 <= 1.0;
subject to dmlaxlon:
	0 <= plaxlon0 + plaxlon6 + plaxlon7 + 13.0*-1 <= 2.0;
subject to dmlaxpar:
	0 <= plaxpar0 + plaxpar6 + plaxpar7 + 8.0*-1 <= 8.0;
subject to dmburlon:
	0 <= pburlon1 + -1 <= 1.0;
subject to dmburpar:
	0 <= pburpar1 + -1 <= 1.0;
subject to dmlonont:
	0 <= plonont1 + -1 <= 1.0;
subject to dmlonoak:
	0 <= plonoak1 + -1 <= 1.0;
subject to dmoakpar:
	0 <= poakpar1 + 2.0*-1 <= 2.0;
subject to dmoaktyo:
	0 <= poaktyo1 + 7.0*-1 <= 7.0;
subject to dmoaktpe:
	0 <= poaktpe1 + 2.0*-1 <= 2.0;
subject to dmontpar:
	0 <= pontpar1 + 2.0*-1 <= 2.0;
subject to dmonttyo:
	0 <= ponttyo1 + 2.0*-1 <= 2.0;
subject to dmonttpe:
	0 <= ponttpe1 + -1 <= 1.0;
subject to dmparsfo:
	0 <= pparsfo1 + 2.0*-1 <= 2.0;
subject to dmrnotyo:
	0 <= prnotyo1 + 5.0*-1 <= 5.0;
subject to dmrnotpe:
	0 <= prnotpe1 + 2.0*-1 <= 2.0;
subject to dmtpeywg:
	0 <= ptpeywg1 + 21.0*-1 <= 3.0;
subject to dmtpeyyz:
	0 <= ptpeyyz1 + 13.0*-1 <= 2.0;
subject to dmtpeyul:
	0 <= ptpeyul1 + 6.0*-1 <= 6.0;
subject to dmtyoyul:
	0 <= ptyoyul1 + 7.0*-1 <= 7.0;
subject to dmtyoyyz:
	0 <= ptyoyyz1 + 17.0*-1 <= 2.0;
subject to dmtyoywg:
	0 <= ptyoywg1 + 25.0*-1 <= 3.0;
subject to dmsfooak:
	-pbosoak6 = 0;
subject to dmlaxbur:
	-pbosbur2 - pburhnl6 = 0;
subject to dmlaxont:
	-pbosont2 - phnlont6 + plaxont0 = 0;
subject to msboshnl:
	0 <= n1007ac1 + n1007ac2 + n1007ac3 + n1047ac1 + n1047ac2 + n1047ac3 + -1;
subject to msbossea:
	0 <= n1007ac1 + n1007ac2 + n1007ac3 + n1032ac1 + n1032ac2 + n1032ac3 + n1032ac4 
	+ n1032ac5 + n1043ac1 + n1043ac2 + n1043ac3 + -1;
subject to msbossfo:
	0 <= n1044ac1 + n1044ac2 + n1044ac3 + n1047ac1 + n1047ac2 + n1047ac3 + 2.0*-1;
subject to msbosyul:
	0 <= n1032ac1 + n1032ac2 + n1032ac3 + n1032ac4 + n1032ac5 + n1033ac1 + n1033ac2 
	+ n1033ac3 + n1033ac4 + n1033ac5 + n1036ac1 + n1036ac2 + n1036ac3 + n1040ac4 + 
	n1040ac5 + n1040ac6 + n1042ac4 + n1042ac5 + n1042ac6 + 7.0*-1;
subject to msbosyvr:
	0 <= n1032ac1 + n1032ac2 + n1032ac3 + n1032ac4 + n1032ac5 + n1033ac1 + n1033ac2 
	+ n1033ac3 + n1033ac4 + n1033ac5 + n1036ac1 + n1036ac2 + n1036ac3 + -1;
subject to msbosywg:
	0 <= n1032ac1 + n1032ac2 + n1032ac3 + n1032ac4 + n1032ac5 + n1033ac1 + n1033ac2 
	+ n1033ac3 + n1033ac4 + n1033ac5 + -1;
subject to msbosyyz:
	0 <= n1032ac1 + n1032ac2 + n1032ac3 + n1032ac4 + n1032ac5 + n1040ac4 + n1040ac5 
	+ n1040ac6 + 3.0*-1;
subject to msburoak:
	0 <= n1005ac3 + n1014ac3 + n1014ac4 + n1014ac5 + n1014ac6 + -1;
subject to msbursea:
	0 <= n1014ac3 + n1014ac4 + n1014ac5 + n1014ac6 + n1016ac3 + n1016ac4 + n1016ac5 
	+ n1016ac6 + -1;
subject to msbursfo:
	0 <= n1016ac3 + n1016ac4 + n1016ac5 + n1016ac6 + 4.0*-1;
subject to mshnllax:
	0 <= n1003ac1 + n1003ac2 + n1003ac3 + n1004ac1 + n1004ac2 + n1004ac3 + n1105ac3 
	+ n1047ac1 + n1047ac2 + n1047ac3 + 3.0*-1;
subject to mshnlsea:
	0 <= n1001ac1 + n1001ac2 + n1001ac3 + n1002ac1 + n1002ac2 + n1002ac3 + n1007ac1 
	+ n1007ac2 + n1007ac3 + 2.0*-1;
subject to mshnlsfo:
	0 <= n1004ac1 + n1004ac2 + n1004ac3 + n1047ac1 + n1047ac2 + n1047ac3 + -1;
subject to mshnlyvr:
	0 <= n1002ac1 + n1002ac2 + n1002ac3 + n1022ac1 + 2.0*-1;
subject to mslassea:
	0 <= n1019ac1 + n1019ac2 + n1019ac3 + n1019ac4 + n1019ac5 + n1021ac1 + n1021ac2 
	+ n1021ac3 + n1021ac4 + n1021ac5 + 4.0*-1;
subject to mslasyvr:
	0 <= n1019ac1 + n1019ac2 + n1019ac3 + n1019ac4 + n1019ac5 + -1;
subject to mslaxoak:
	0 <= n1105ac3 + n1017ac3 + n1017ac4 + n1017ac5 + n1017ac6 + n1046ac3 + 2.0*-1;
subject to mslaxsea:
	0 <= n1008ac1 + n1008ac2 + n1008ac3 + n1008ac4 + n1008ac5 + n1008ac6 + n1009ac1 
	+ n1009ac2 + n1009ac3 + n1009ac4 + n1009ac5 + n1011ac1 + n1011ac2 + n1011ac3 + 
	n1011ac4 + n1011ac5 + n1011ac6 + n1012ac1 + n1012ac2 + n1012ac3 + n1012ac4 + 
	n1012ac5 + n1017ac3 + n1017ac4 + n1017ac5 + n1017ac6 + n1026ac1 + n1027ac1 + 
	n1029ac1 + n1030ac1 + n1050ac3 + n1050ac4 + n1050ac5 + 7.0*-1;
subject to mslaxsfo:
	0 <= n1004ac1 + n1004ac2 + n1004ac3 + n1008ac1 + n1008ac2 + n1008ac3 + n1008ac4 
	+ n1008ac5 + n1008ac6 + n1011ac1 + n1011ac2 + n1011ac3 + n1011ac4 + n1011ac5 + 
	n1011ac6 + n1029ac1 + n1030ac1 + n1047ac1 + n1047ac2 + n1047ac3 + n1051ac1 + 
	n1051ac2 + n1051ac3 + n1051ac4 + n1051ac5 + n1051ac6 + 21.0*-1;
subject to mslaxtpe:
	0 <= n1027ac1 + n1029ac1 + n1030ac1 + 2.0*-1;
subject to mslaxyvr:
	0 <= n1008ac1 + n1008ac2 + n1008ac3 + n1008ac4 + n1008ac5 + n1008ac6 + n1009ac1 
	+ n1009ac2 + n1009ac3 + n1009ac4 + n1009ac5 + n1029ac1 + 3.0*-1;
subject to mslonpar:
	0 <= n1022ac1 + n1023ac1 + n1026ac1 + -1;
subject to mslonsea:
	0 <= n1023ac1 + n1026ac1 + -1;
subject to mslonyvr:
	0 <= n1022ac1 + n1023ac1 + -1;
subject to msoakont:
	0 <= n1006ac3 + n1013ac3 + n1013ac4 + n1013ac5 + n1013ac6 + -1;
subject to msoaksea:
	0 <= n1013ac3 + n1013ac4 + n1013ac5 + n1013ac6 + n1014ac3 + n1014ac4 + n1014ac5 
	+ n1014ac6 + n1017ac3 + n1017ac4 + n1017ac5 + n1017ac6 + 3.0*-1;
subject to msontsfo:
	0 <= n1015ac3 + n1015ac4 + n1015ac5 + n1015ac6 + 3.0*-1;
subject to msontsea:
	0 <= n1013ac3 + n1013ac4 + n1013ac5 + n1013ac6 + n1015ac3 + n1015ac4 + n1015ac5 
	+ n1015ac6 + n1050ac3 + n1050ac4 + n1050ac5 + 2.0*-1;
subject to msparsea:
	0 >= n1023ac1 + n1026ac1 + -1;
subject to msrnosea:
	0 <= n1018ac1 + n1018ac2 + n1018ac3 + n1018ac4 + n1018ac5 + n1018ac6 + n1020ac1 
	+ n1020ac2 + n1020ac3 + n1020ac4 + n1020ac5 + n1020ac6 + 4.0*-1;
subject to msrnoyvr:
	0 <= n1018ac1 + n1018ac2 + n1018ac3 + n1018ac4 + n1018ac5 + n1018ac6 + -1;
subject to msseasfo:
	0 <= n1008ac1 + n1008ac2 + n1008ac3 + n1008ac4 + n1008ac5 + n1008ac6 + n1011ac1 
	+ n1011ac2 + n1011ac3 + n1011ac4 + n1011ac5 + n1011ac6 + n1015ac3 + n1015ac4 + 
	n1015ac5 + n1015ac6 + n1016ac3 + n1016ac4 + n1016ac5 + n1016ac6 + n1028ac1 + 
	n1029ac1 + n1030ac1 + 10.0*-1;
subject to msseatpe:
	0 <= n1027ac1 + n1028ac1 + n1029ac1 + n1030ac1 + -1;
subject to msseatyo:
	0 <= n1027ac1 + n1028ac1 + n1029ac1 + n1030ac1 + -1;
subject to msseayvr:
	0 <= n1002ac1 + n1002ac2 + n1002ac3 + n1008ac1 + n1008ac2 + n1008ac3 + n1008ac4 
	+ n1008ac5 + n1008ac6 + n1009ac1 + n1009ac2 + n1009ac3 + n1009ac4 + n1009ac5 + 
	n1010ac1 + n1010ac2 + n1010ac3 + n1010ac4 + n1010ac5 + n1010ac6 + n1018ac1 + 
	n1018ac2 + n1018ac3 + n1018ac4 + n1018ac5 + n1018ac6 + n1019ac1 + n1019ac2 + 
	n1019ac3 + n1019ac4 + n1019ac5 + n1023ac1 + n1029ac1 + n1032ac1 + n1032ac2 + 
	n1032ac3 + n1032ac4 + n1032ac5 + 6.0*-1;
subject to mssfoyvr:
	0 <= n1008ac1 + n1008ac2 + n1008ac3 + n1008ac4 + n1008ac5 + n1008ac6 + n1029ac1 
	+ 5.0*-1;
subject to mstpetyo:
	0 <= n1027ac1 + n1028ac1 + n1029ac1 + n1030ac1 + -1;
subject to msyulyvr:
	0 <= n1032ac1 + n1032ac2 + n1032ac3 + n1032ac4 + n1032ac5 + n1033ac1 + n1033ac2 
	+ n1033ac3 + n1033ac4 + n1033ac5 + n1034ac1 + n1034ac2 + n1034ac3 + n1036ac1 + 
	n1036ac2 + n1036ac3 + 7.0*-1;
subject to msyulywg:
	0 <= n1032ac1 + n1032ac2 + n1032ac3 + n1032ac4 + n1032ac5 + n1033ac1 + n1033ac2 
	+ n1033ac3 + n1033ac4 + n1033ac5 + n1037ac4 + n1037ac5 + n1039ac4 + n1039ac5 + 
	5.0*-1;
subject to msyulyyz:
	0 <= n1032ac1 + n1032ac2 + n1032ac3 + n1032ac4 + n1032ac5 + n1034ac1 + n1034ac2 
	+ n1034ac3 + n1037ac4 + n1037ac5 + n1040ac4 + n1040ac5 + n1040ac6 + n1041ac4 + 
	n1041ac5 + n1041ac6 + 24.0*-1;
subject to msyvrywg:
	0 <= n1032ac1 + n1032ac2 + n1032ac3 + n1032ac4 + n1032ac5 + n1033ac1 + n1033ac2 
	+ n1033ac3 + n1033ac4 + n1033ac5 + n1035ac1 + n1035ac2 + n1035ac3 + n1035ac4 + 
	n1035ac5 + 5.0*-1;
subject to msyvryyz:
	0 <= n1032ac1 + n1032ac2 + n1032ac3 + n1032ac4 + n1032ac5 + n1034ac1 + n1034ac2 
	+ n1034ac3 + n1035ac1 + n1035ac2 + n1035ac3 + n1035ac4 + n1035ac5 + 5.0*-1;
subject to msywgyyz:
	0 <= n1032ac1 + n1032ac2 + n1032ac3 + n1032ac4 + n1032ac5 + n1035ac1 + n1035ac2 
	+ n1035ac3 + n1035ac4 + n1035ac5 + n1037ac4 + n1037ac5 + n1038ac4 + n1038ac5 + 
	4.0*-1;
subject to msbosoak:
	0 <= n1046ac3 + -1;
subject to mshnloak:
	0 <= n1005ac3 + n1105ac3 + n1006ac3 + -1;
subject to mslaxlon:
	0 <= n1026ac1 + -1;
subject to p1044x32:
	0 <= -pbossfo0 + 73.0*n1044ac1 + 73.0*n1044ac2 + 73.0*n1044ac3;
subject to p1047x54:
	0 <= -pbossfo1 + 73.0*n1047ac1 + 73.0*n1047ac2 + 73.0*n1047ac3;
subject to p1032x76:
	0 <= -pbosyul0 + 115.0*n1032ac1 + 115.0*n1032ac2 + 115.0*n1032ac3 + 
	115.0*n1032ac4 + 115.0*n1032ac5;
subject to p1033x54:
	0 <= -pbosyul1 + 115.0*n1033ac1 + 115.0*n1033ac2 + 115.0*n1033ac3 + 
	115.0*n1033ac4 + 115.0*n1033ac5;
subject to p1036x43:
	0 <= -pbosyul2 + 115.0*n1036ac1 + 115.0*n1036ac2 + 115.0*n1036ac3;
subject to p1040x43:
	0 <= -pbosyul3 + 115.0*n1040ac4 + 115.0*n1040ac5 + 115.0*n1040ac6;
subject to p1042x32:
	0 <= -pbosyul4 + 115.0*n1042ac4 + 115.0*n1042ac5 + 115.0*n1042ac6;
subject to p1032x75:
	0 <= -pbosyyz0 + 85.0*n1032ac1 + 85.0*n1032ac2 + 85.0*n1032ac3 + 85.0*n1032ac4 
	+ 85.0*n1032ac5;
subject to p1040x42:
	0 <= -pbosyyz1 + 85.0*n1040ac4 + 85.0*n1040ac5 + 85.0*n1040ac6;
subject to p1016x43:
	0 <= -pbursfo0 + 81.0*n1016ac3 + 81.0*n1016ac4 + 81.0*n1016ac5 + 81.0*n1016ac6;
subject to p1003x32:
	0 <= -phnllax0 + 118.0*n1003ac1 + 118.0*n1003ac2 + 118.0*n1003ac3;
subject to p1004x43:
	0 <= -phnllax1 + 118.0*n1004ac1 + 118.0*n1004ac2 + 118.0*n1004ac3;
subject to p1105x43:
	0 <= -phnllax2 + 118.0*n1105ac3;
subject to p1047x23:
	0 <= -phnllax3 + 118.0*n1047ac1 + 118.0*n1047ac2 + 118.0*n1047ac3;
subject to p1001x32:
	0 <= -phnlsea0 + 67.0*n1001ac1 + 67.0*n1001ac2 + 67.0*n1001ac3;
subject to p1002x43:
	0 <= -phnlsea1 + 67.0*n1002ac1 + 67.0*n1002ac2 + 67.0*n1002ac3;
subject to p1007x43:
	0 <= -phnlsea2 + 67.0*n1007ac1 + 67.0*n1007ac2 + 67.0*n1007ac3;
subject to p1002x42:
	0 <= -phnlyvr0 + 40.0*n1002ac1 + 40.0*n1002ac2 + 40.0*n1002ac3;
subject to p1022x23:
	0 <= -phnlyvr1 + 40.0*n1022ac1;
subject to p1019x43:
	0 <= -plassea0 + 110.0*n1019ac1 + 110.0*n1019ac2 + 110.0*n1019ac3 + 
	110.0*n1019ac4 + 110.0*n1019ac5;
subject to p1021x32:
	0 <= -plassea1 + 110.0*n1021ac1 + 110.0*n1021ac2 + 110.0*n1021ac3 + 
	110.0*n1021ac4 + 110.0*n1021ac5;
subject to p1105x32:
	0 <= -plaxoak0 + 46.0*n1105ac3;
subject to p1017x43:
	0 <= -plaxoak1 + 46.0*n1017ac3 + 46.0*n1017ac4 + 46.0*n1017ac5 + 46.0*n1017ac6;
subject to p1046x23:
	0 <= -plaxoak2 + 46.0*n1046ac3;
subject to p1008x53:
	0 <= -plaxsea0 + 139.0*n1008ac1 + 139.0*n1008ac2 + 139.0*n1008ac3 + 
	139.0*n1008ac4 + 139.0*n1008ac5 + 139.0*n1008ac6;
subject to p1009x43:
	0 <= -plaxsea1 + 139.0*n1009ac1 + 139.0*n1009ac2 + 139.0*n1009ac3 + 
	139.0*n1009ac4 + 139.0*n1009ac5;
subject to p1011x42:
	0 <= -plaxsea2 + 139.0*n1011ac1 + 139.0*n1011ac2 + 139.0*n1011ac3 + 
	139.0*n1011ac4 + 139.0*n1011ac5 + 139.0*n1011ac6;
subject to p1012x32:
	0 <= -plaxsea3 + 139.0*n1012ac1 + 139.0*n1012ac2 + 139.0*n1012ac3 + 
	139.0*n1012ac4 + 139.0*n1012ac5;
subject to p1017x42:
	0 <= -plaxsea4 + 139.0*n1017ac3 + 139.0*n1017ac4 + 139.0*n1017ac5 + 
	139.0*n1017ac6;
subject to p1026x23:
	0 <= -plaxsea5 + 139.0*n1026ac1;
subject to p1027x23:
	0 <= -plaxsea6 + 139.0*n1027ac1;
subject to p1029x24:
	0 <= -plaxsea7 + 139.0*n1029ac1;
subject to p1030x24:
	0 <= -plaxsea8 + 139.0*n1030ac1;
subject to p1050x32:
	0 <= -plaxsea9 + 139.0*n1050ac3 + 139.0*n1050ac4 + 139.0*n1050ac5;
subject to p1004x32:
	0 <= -plaxsfo0 + 168.0*n1004ac1 + 168.0*n1004ac2 + 168.0*n1004ac3;
subject to p1008x54:
	0 <= -plaxsfo1 + 168.0*n1008ac1 + 168.0*n1008ac2 + 168.0*n1008ac3 + 
	168.0*n1008ac4 + 168.0*n1008ac5 + 168.0*n1008ac6;
subject to p1011x43:
	0 <= -plaxsfo2 + 168.0*n1011ac1 + 168.0*n1011ac2 + 168.0*n1011ac3 + 
	168.0*n1011ac4 + 168.0*n1011ac5 + 168.0*n1011ac6;
subject to p1029x23:
	0 <= -plaxsfo3 + 168.0*n1029ac1;
subject to p1030x23:
	0 <= -plaxsfo4 + 168.0*n1030ac1;
subject to p1047x34:
	0 <= -plaxsfo5 + 168.0*n1047ac1 + 168.0*n1047ac2 + 168.0*n1047ac3;
subject to p1051x23:
	0 <= -plaxsfo6 + 168.0*n1051ac1 + 168.0*n1051ac2 + 168.0*n1051ac3 + 
	168.0*n1051ac4 + 168.0*n1051ac5 + 168.0*n1051ac6;
subject to p1027x25:
	0 <= -plaxtpe0 + 18.0*n1027ac1;
subject to p1029x27:
	0 <= -plaxtpe1 + 18.0*n1029ac1;
subject to p1030x26:
	0 <= -plaxtpe2 + 18.0*n1030ac1;
subject to p1008x52:
	0 <= -plaxyvr0 + 77.0*n1008ac1 + 77.0*n1008ac2 + 77.0*n1008ac3 + 77.0*n1008ac4 
	+ 77.0*n1008ac5 + 77.0*n1008ac6;
subject to p1009x42:
	0 <= -plaxyvr1 + 77.0*n1009ac1 + 77.0*n1009ac2 + 77.0*n1009ac3 + 77.0*n1009ac4 
	+ 77.0*n1009ac5;
subject to p1029x25:
	0 <= -plaxyvr2 + 77.0*n1029ac1;
subject to p1013x32:
	0 <= -poaksea0 + 43.0*n1013ac3 + 43.0*n1013ac4 + 43.0*n1013ac5 + 43.0*n1013ac6;
subject to p1014x32:
	0 <= -poaksea1 + 43.0*n1014ac3 + 43.0*n1014ac4 + 43.0*n1014ac5 + 43.0*n1014ac6;
subject to p1017x32:
	0 <= -poaksea2 + 43.0*n1017ac3 + 43.0*n1017ac4 + 43.0*n1017ac5 + 43.0*n1017ac6;
subject to p1015x43:
	0 <= -pontsfo0 + 69.0*n1015ac3 + 69.0*n1015ac4 + 69.0*n1015ac5 + 69.0*n1015ac6;
subject to p1013x42:
	0 <= -pontsea0 + 25.0*n1013ac3 + 25.0*n1013ac4 + 25.0*n1013ac5 + 25.0*n1013ac6;
subject to p1015x42:
	0 <= -pontsea1 + 25.0*n1015ac3 + 25.0*n1015ac4 + 25.0*n1015ac5 + 25.0*n1015ac6;
subject to p1050x42:
	0 <= -pontsea2 + 25.0*n1050ac3 + 25.0*n1050ac4 + 25.0*n1050ac5;
subject to p1018x43:
	0 <= -prnosea0 + 85.0*n1018ac1 + 85.0*n1018ac2 + 85.0*n1018ac3 + 85.0*n1018ac4 
	+ 85.0*n1018ac5 + 85.0*n1018ac6;
subject to p1020x32:
	0 <= -prnosea1 + 85.0*n1020ac1 + 85.0*n1020ac2 + 85.0*n1020ac3 + 85.0*n1020ac4 
	+ 85.0*n1020ac5 + 85.0*n1020ac6;
subject to p1008x34:
	0 <= -pseasfo0 + 170.0*n1008ac1 + 170.0*n1008ac2 + 170.0*n1008ac3 + 
	170.0*n1008ac4 + 170.0*n1008ac5 + 170.0*n1008ac6;
subject to p1011x23:
	0 <= -pseasfo1 + 170.0*n1011ac1 + 170.0*n1011ac2 + 170.0*n1011ac3 + 
	170.0*n1011ac4 + 170.0*n1011ac5 + 170.0*n1011ac6;
subject to p1015x23:
	0 <= -pseasfo2 + 170.0*n1015ac3 + 170.0*n1015ac4 + 170.0*n1015ac5 + 
	170.0*n1015ac6;
subject to p1016x23:
	0 <= -pseasfo3 + 170.0*n1016ac3 + 170.0*n1016ac4 + 170.0*n1016ac5 + 
	170.0*n1016ac6;
subject to p1028x32:
	0 <= -pseasfo4 + 170.0*n1028ac1;
subject to p1029x43:
	0 <= -pseasfo5 + 170.0*n1029ac1;
subject to p1030x43:
	0 <= -pseasfo6 + 170.0*n1030ac1;
subject to p1002x32:
	0 <= -pseayvr0 + 109.0*n1002ac1 + 109.0*n1002ac2 + 109.0*n1002ac3;
subject to p1008x32:
	0 <= -pseayvr1 + 109.0*n1008ac1 + 109.0*n1008ac2 + 109.0*n1008ac3 + 
	109.0*n1008ac4 + 109.0*n1008ac5 + 109.0*n1008ac6;
subject to p1009x32:
	0 <= -pseayvr2 + 109.0*n1009ac1 + 109.0*n1009ac2 + 109.0*n1009ac3 + 
	109.0*n1009ac4 + 109.0*n1009ac5;
subject to p1010x32:
	0 <= -pseayvr3 + 109.0*n1010ac1 + 109.0*n1010ac2 + 109.0*n1010ac3 + 
	109.0*n1010ac4 + 109.0*n1010ac5 + 109.0*n1010ac6;
subject to p1018x32:
	0 <= -pseayvr4 + 109.0*n1018ac1 + 109.0*n1018ac2 + 109.0*n1018ac3 + 
	109.0*n1018ac4 + 109.0*n1018ac5 + 109.0*n1018ac6;
subject to p1019x32:
	0 <= -pseayvr5 + 109.0*n1019ac1 + 109.0*n1019ac2 + 109.0*n1019ac3 + 
	109.0*n1019ac4 + 109.0*n1019ac5;
subject to p1023x32:
	0 <= -pseayvr6 + 109.0*n1023ac1;
subject to p1029x45:
	0 <= -pseayvr7 + 109.0*n1029ac1;
subject to p1032x23:
	0 <= -pseayvr8 + 109.0*n1032ac1 + 109.0*n1032ac2 + 109.0*n1032ac3 + 
	109.0*n1032ac4 + 109.0*n1032ac5;
subject to p1008x42:
	0 <= -psfoyvr0 + 71.0*n1008ac1 + 71.0*n1008ac2 + 71.0*n1008ac3 + 71.0*n1008ac4 
	+ 71.0*n1008ac5 + 71.0*n1008ac6;
subject to p1029x35:
	0 <= -psfoyvr1 + 71.0*n1029ac1;
subject to p1032x63:
	0 <= -pyulyvr0 + 44.0*n1032ac1 + 44.0*n1032ac2 + 44.0*n1032ac3 + 44.0*n1032ac4 
	+ 44.0*n1032ac5;
subject to p1033x42:
	0 <= -pyulyvr1 + 44.0*n1033ac1 + 44.0*n1033ac2 + 44.0*n1033ac3 + 44.0*n1033ac4 
	+ 44.0*n1033ac5;
subject to p1034x42:
	0 <= -pyulyvr2 + 44.0*n1034ac1 + 44.0*n1034ac2 + 44.0*n1034ac3;
subject to p1036x32:
	0 <= -pyulyvr3 + 44.0*n1036ac1 + 44.0*n1036ac2 + 44.0*n1036ac3;
subject to p1032x64:
	0 <= -pyulywg0 + 99.0*n1032ac1 + 99.0*n1032ac2 + 99.0*n1032ac3 + 99.0*n1032ac4 
	+ 99.0*n1032ac5;
subject to p1033x43:
	0 <= -pyulywg1 + 99.0*n1033ac1 + 99.0*n1033ac2 + 99.0*n1033ac3 + 99.0*n1033ac4 
	+ 99.0*n1033ac5;
subject to p1037x42:
	0 <= -pyulywg2 + 99.0*n1037ac4 + 99.0*n1037ac5;
subject to p1039x32:
	0 <= -pyulywg3 + 99.0*n1039ac4 + 99.0*n1039ac5;
subject to p1032x65:
	0 <= -pyulyyz0 + 130.0*n1032ac1 + 130.0*n1032ac2 + 130.0*n1032ac3 + 
	130.0*n1032ac4 + 130.0*n1032ac5;
subject to p1034x43:
	0 <= -pyulyyz1 + 130.0*n1034ac1 + 130.0*n1034ac2 + 130.0*n1034ac3;
subject to p1037x43:
	0 <= -pyulyyz2 + 130.0*n1037ac4 + 130.0*n1037ac5;
subject to p1040x32:
	0 <= -pyulyyz3 + 130.0*n1040ac4 + 130.0*n1040ac5 + 130.0*n1040ac6;
subject to p1041x32:
	0 <= -pyulyyz4 + 130.0*n1041ac4 + 130.0*n1041ac5 + 130.0*n1041ac6;
subject to p1032x34:
	0 <= -pyvrywg0 + 89.0*n1032ac1 + 89.0*n1032ac2 + 89.0*n1032ac3 + 89.0*n1032ac4 
	+ 89.0*n1032ac5;
subject to p1033x23:
	0 <= -pyvrywg1 + 89.0*n1033ac1 + 89.0*n1033ac2 + 89.0*n1033ac3 + 89.0*n1033ac4 
	+ 89.0*n1033ac5;
subject to p1035x23:
	0 <= -pyvrywg2 + 89.0*n1035ac1 + 89.0*n1035ac2 + 89.0*n1035ac3 + 89.0*n1035ac4 
	+ 89.0*n1035ac5;
subject to p1032x35:
	0 <= -pyvryyz0 + 76.0*n1032ac1 + 76.0*n1032ac2 + 76.0*n1032ac3 + 76.0*n1032ac4 
	+ 76.0*n1032ac5;
subject to p1034x23:
	0 <= -pyvryyz1 + 76.0*n1034ac1 + 76.0*n1034ac2 + 76.0*n1034ac3;
subject to p1035x24:
	0 <= -pyvryyz2 + 76.0*n1035ac1 + 76.0*n1035ac2 + 76.0*n1035ac3 + 76.0*n1035ac4 
	+ 76.0*n1035ac5;
subject to p1032x45:
	0 <= -pywgyyz0 + 83.0*n1032ac1 + 83.0*n1032ac2 + 83.0*n1032ac3 + 83.0*n1032ac4 
	+ 83.0*n1032ac5;
subject to p1035x34:
	0 <= -pywgyyz1 + 83.0*n1035ac1 + 83.0*n1035ac2 + 83.0*n1035ac3 + 83.0*n1035ac4 
	+ 83.0*n1035ac5;
subject to p1037x23:
	0 <= -pywgyyz2 + 83.0*n1037ac4 + 83.0*n1037ac5;
subject to p1038x23:
	0 <= -pywgyyz3 + 83.0*n1038ac4 + 83.0*n1038ac5;

solve;
	display pboshnl0;
	display pboshnl1;
	display pboshnl7;
	display pboshnl8;
	display pboslax0;
	display pboslax1;
	display pboslax7;
	display pbossea0;
	display pbossea1;
	display pbossea2;
	display pbossfo0;
	display pbossfo1;
	display pbostpe1;
	display pbostpe2;
	display pbostyo1;
	display pbostyo2;
	display pbosyul0;
	display pbosyul1;
	display pbosyul2;
	display pbosyul3;
	display pbosyul4;
	display pbosyvr0;
	display pbosyvr1;
	display pbosyvr2;
	display pbosywg0;
	display pbosywg1;
	display pbosywg7;
	display pbosyyz0;
	display pbosyyz1;
	display pburoak0;
	display pburoak1;
	display pbursea0;
	display pbursea1;
	display pbursfo0;
	display phnllax0;
	display phnllax1;
	display phnllax2;
	display phnllax3;
	display phnllon0;
	display phnllon6;
	display phnlpar0;
	display phnlpar6;
	display phnlsea0;
	display phnlsea1;
	display phnlsea2;
	display phnlsfo0;
	display phnlsfo1;
	display phnlsfo7;
	display phnlyvr0;
	display phnlyvr1;
	display phnlyvr7;
	display plassea0;
	display plassea1;
	display plasyvr0;
	display plasyvr6;
	display plaxoak0;
	display plaxoak1;
	display plaxoak2;
	display plaxsea0;
	display plaxsea1;
	display plaxsea2;
	display plaxsea3;
	display plaxsea4;
	display plaxsea5;
	display plaxsea6;
	display plaxsea7;
	display plaxsea8;
	display plaxsea9;
	display plaxsfo0;
	display plaxsfo1;
	display plaxsfo2;
	display plaxsfo3;
	display plaxsfo4;
	display plaxsfo5;
	display plaxsfo6;
	display plaxtpe0;
	display plaxtpe1;
	display plaxtpe2;
	display plaxtpe8;
	display plaxtyo0;
	display plaxtyo1;
	display plaxtyo2;
	display plaxtyo8;
	display plaxyvr0;
	display plaxyvr1;
	display plaxyvr2;
	display plonpar0;
	display plonpar1;
	display plonpar2;
	display plonsea0;
	display plonsea1;
	display plonyvr0;
	display plonyvr1;
	display plonyvr7;
	display poakont0;
	display poakont1;
	display poaksea0;
	display poaksea1;
	display poaksea2;
	display pontsfo0;
	display pontsea0;
	display pontsea1;
	display pontsea2;
	display pparsea0;
	display pparsea1;
	display pparyvr0;
	display pparyvr1;
	display pparyvr7;
	display prnosea0;
	display prnosea1;
	display prnoyvr0;
	display pseasfo0;
	display pseasfo1;
	display pseasfo2;
	display pseasfo3;
	display pseasfo4;
	display pseasfo5;
	display pseasfo6;
	display pseatpe0;
	display pseatpe1;
	display pseatpe2;
	display pseatpe3;
	display pseatyo0;
	display pseatyo1;
	display pseatyo2;
	display pseatyo3;
	display pseayvr0;
	display pseayvr1;
	display pseayvr2;
	display pseayvr3;
	display pseayvr4;
	display pseayvr5;
	display pseayvr6;
	display pseayvr7;
	display pseayvr8;
	display psfotpe0;
	display psfotpe1;
	display psfotpe2;
	display psfotpe8;
	display psfotyo0;
	display psfotyo1;
	display psfotyo2;
	display psfotyo8;
	display psfoyvr0;
	display psfoyvr1;
	display ptpetyo0;
	display ptpetyo1;
	display ptpetyo2;
	display ptpetyo3;
	display ptpeyvr0;
	display ptyoyvr0;
	display pyulyvr0;
	display pyulyvr1;
	display pyulyvr2;
	display pyulyvr3;
	display pyulywg0;
	display pyulywg1;
	display pyulywg2;
	display pyulywg3;
	display pyulyyz0;
	display pyulyyz1;
	display pyulyyz2;
	display pyulyyz3;
	display pyulyyz4;
	display pyvrywg0;
	display pyvrywg1;
	display pyvrywg2;
	display pyvryyz0;
	display pyvryyz1;
	display pyvryyz2;
	display pywgyyz0;
	display pywgyyz1;
	display pywgyyz2;
	display pywgyyz3;
	display pbosoak0;
	display pbosoak6;
	display pbosbur1;
	display pbosbur2;
	display pbosont1;
	display pbosont2;
	display pburyvr1;
	display pburtyo1;
	display pburtpe1;
	display pburhnl0;
	display pburhnl6;
	display phnloak0;
	display phnloak1;
	display phnloak2;
	display phnloak8;
	display phnlont0;
	display phnlont6;
	display phnlywg1;
	display phnlyyz1;
	display phnlyul1;
	display plastyo1;
	display plastpe1;
	display plaxlon0;
	display plaxlon6;
	display plaxlon7;
	display plaxpar0;
	display plaxpar6;
	display plaxpar7;
	display pburlon1;
	display pburpar1;
	display plonont1;
	display plonoak1;
	display poakpar1;
	display poaktyo1;
	display poaktpe1;
	display pontpar1;
	display ponttyo1;
	display ponttpe1;
	display pparsfo1;
	display prnotyo1;
	display prnotpe1;
	display ptpeywg1;
	display ptpeyyz1;
	display ptpeyul1;
	display ptyoyul1;
	display ptyoyyz1;
	display ptyoywg1;
	display plaxont0;
	display grdtimo1;
	display grdtimn1;
	display grdtimo2;
	display grdtimn2;
	display grdtimo3;
	display grdtimn3;
	display grdtimo4;
	display grdtimn4;
	display grdtimo5;
	display grdtimn5;
	display grdtimo6;
	display grdtimn6;
	display n1001ac1;
	display n1001ac2;
	display n1001ac3;
	display n1002ac1;
	display n1002ac2;
	display n1002ac3;
	display n1003ac1;
	display n1003ac2;
	display n1003ac3;
	display n1004ac1;
	display n1004ac2;
	display n1004ac3;
	display n1005ac3;
	display n1105ac3;
	display n1006ac3;
	display n1007ac1;
	display n1007ac2;
	display n1007ac3;
	display n1008ac1;
	display n1008ac2;
	display n1008ac3;
	display n1008ac4;
	display n1008ac5;
	display n1008ac6;
	display n1009ac1;
	display n1009ac2;
	display n1009ac3;
	display n1009ac4;
	display n1009ac5;
	display n1010ac1;
	display n1010ac2;
	display n1010ac3;
	display n1010ac4;
	display n1010ac5;
	display n1010ac6;
	display n1011ac1;
	display n1011ac2;
	display n1011ac3;
	display n1011ac4;
	display n1011ac5;
	display n1011ac6;
	display n1012ac1;
	display n1012ac2;
	display n1012ac3;
	display n1012ac4;
	display n1012ac5;
	display n1013ac3;
	display n1013ac4;
	display n1013ac5;
	display n1013ac6;
	display n1014ac3;
	display n1014ac4;
	display n1014ac5;
	display n1014ac6;
	display n1015ac3;
	display n1015ac4;
	display n1015ac5;
	display n1015ac6;
	display n1016ac3;
	display n1016ac4;
	display n1016ac5;
	display n1016ac6;
	display n1017ac3;
	display n1017ac4;
	display n1017ac5;
	display n1017ac6;
	display n1018ac1;
	display n1018ac2;
	display n1018ac3;
	display n1018ac4;
	display n1018ac5;
	display n1018ac6;
	display n1019ac1;
	display n1019ac2;
	display n1019ac3;
	display n1019ac4;
	display n1019ac5;
	display n1020ac1;
	display n1020ac2;
	display n1020ac3;
	display n1020ac4;
	display n1020ac5;
	display n1020ac6;
	display n1021ac1;
	display n1021ac2;
	display n1021ac3;
	display n1021ac4;
	display n1021ac5;
	display n1022ac1;
	display n1023ac1;
	display n1026ac1;
	display n1027ac1;
	display n1028ac1;
	display n1029ac1;
	display n1030ac1;
	display n1032ac1;
	display n1032ac2;
	display n1032ac3;
	display n1032ac4;
	display n1032ac5;
	display n1033ac1;
	display n1033ac2;
	display n1033ac3;
	display n1033ac4;
	display n1033ac5;
	display n1034ac1;
	display n1034ac2;
	display n1034ac3;
	display n1035ac1;
	display n1035ac2;
	display n1035ac3;
	display n1035ac4;
	display n1035ac5;
	display n1036ac1;
	display n1036ac2;
	display n1036ac3;
	display n1037ac4;
	display n1037ac5;
	display n1038ac4;
	display n1038ac5;
	display n1039ac4;
	display n1039ac5;
	display n1040ac4;
	display n1040ac5;
	display n1040ac6;
	display n1041ac4;
	display n1041ac5;
	display n1041ac6;
	display n1042ac4;
	display n1042ac5;
	display n1042ac6;
	display n1043ac1;
	display n1043ac2;
	display n1043ac3;
	display n1044ac1;
	display n1044ac2;
	display n1044ac3;
	display n1046ac3;
	display n1047ac1;
	display n1047ac2;
	display n1047ac3;
	display n1050ac3;
	display n1050ac4;
	display n1050ac5;
	display n1051ac1;
	display n1051ac2;
	display n1051ac3;
	display n1051ac4;
	display n1051ac5;
	display n1051ac6;
display obj;
