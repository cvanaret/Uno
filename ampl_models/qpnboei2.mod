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

#   Source: a variant on the BOEING2 linear programming problem
#   with an additional NONCONVEX diagonal Hessian matrix as given by
#   N. I. M. Gould, "An algorithm for large-scale quadratic programming",
#   IMA J. Num. Anal (1991), 11, 299-324, problem class 4.

#   SIF input: Nick Gould, January 1993

#   classification QLR2-MN-143-166

#
#   Problem :
#   *********
#   Source: a variant on the BOEING2 linear programming problem
#   with an additional NONCONVEX diagonal Hessian matrix as given by
#   N. I. M. Gould, "An algorithm for large-scale quadratic programming",
#   IMA J. Num. Anal (1991), 11, 299-324, problem class 4.
#   SIF input: Nick Gould, January 1993
	param n := 143;

	var pbosord0 >= 0.0, := 0.0;
	var pbosord1 >= 0.0, := 0.0;
	var pbosord2 >= 0.0, := 0.0;
	var pbosord3 >= 0.0, := 0.0;
	var pbosord4 >= 0.0, := 0.0;
	var pboslga0 >= 0.0, := 0.0;
	var pboslga1 >= 0.0, := 0.0;
	var pboslga2 >= 0.0, := 0.0;
	var pboslga3 >= 0.0, := 0.0;
	var pboscle0 >= 0.0, := 0.0;
	var pboscle1 >= 0.0, := 0.0;
	var pboscle2 >= 0.0, := 0.0;
	var pboscle3 >= 0.0, := 0.0;
	var pordbos0 >= 0.0, := 0.0;
	var pordbos1 >= 0.0, := 0.0;
	var pordbos2 >= 0.0, := 0.0;
	var pordbos3 >= 0.0, := 0.0;
	var pordbos4 >= 0.0, := 0.0;
	var pordlga0 >= 0.0, := 0.0;
	var pordlga1 >= 0.0, := 0.0;
	var pordlga2 >= 0.0, := 0.0;
	var pordlga3 >= 0.0, := 0.0;
	var pordcle0 >= 0.0, := 0.0;
	var pordcle1 >= 0.0, := 0.0;
	var pordcle2 >= 0.0, := 0.0;
	var pordcle3 >= 0.0, := 0.0;
	var plgabos0 >= 0.0, := 0.0;
	var plgabos1 >= 0.0, := 0.0;
	var plgabos2 >= 0.0, := 0.0;
	var plgabos3 >= 0.0, := 0.0;
	var plgabos4 >= 0.0, := 0.0;
	var plgabos5 >= 0.0, := 0.0;
	var plgaord0 >= 0.0, := 0.0;
	var plgaord1 >= 0.0, := 0.0;
	var plgaord2 >= 0.0, := 0.0;
	var plgaord3 >= 0.0, := 0.0;
	var plgacle0 >= 0.0, := 0.0;
	var plgacle1 >= 0.0, := 0.0;
	var plgacle2 >= 0.0, := 0.0;
	var pclebos0 >= 0.0, := 0.0;
	var pclebos1 >= 0.0, := 0.0;
	var pclebos2 >= 0.0, := 0.0;
	var pclebos3 >= 0.0, := 0.0;
	var pclebos4 >= 0.0, := 0.0;
	var pcleord0 >= 0.0, := 0.0;
	var pcleord1 >= 0.0, := 0.0;
	var pcleord2 >= 0.0, := 0.0;
	var pcleord3 >= 0.0, := 0.0;
	var pcleord4 >= 0.0, := 0.0;
	var pclelga0 >= 0.0, := 0.0;
	var pclelga1 >= 0.0, := 0.0;
	var pclelga2 >= 0.0, := 0.0;
	var pclelga3 >= 0.0, := 0.0;
	var bbosord0 >= 0.0, := 0.0;
	var bbosord1 >= 0.0, := 0.0;
	var bboscle0 >= 0.0, := 0.0;
	var bboscle1 >= 0.0, := 0.0;
	var bordbos0 >= 0.0, := 0.0;
	var cbosord0 >= 0.0, := 0.0;
	var cbosord1 >= 0.0, := 0.0;
	var cbosord2 >= 0.0, := 0.0;
	var cbosord3 >= 0.0, := 0.0;
	var cbosord4 >= 0.0, := 0.0;
	var cboscle0 >= 0.0, := 0.0;
	var cboscle1 >= 0.0, := 0.0;
	var cboscle2 >= 0.0, := 0.0;
	var cboscle3 >= 0.0, := 0.0;
	var cordbos0 >= 0.0, := 0.0;
	var cordbos1 >= 0.0, := 0.0;
	var cordbos2 >= 0.0, := 0.0;
	var cordbos3 >= 0.0, := 0.0;
	var cordbos4 >= 0.0, := 0.0;
	var cordlga0 >= 0.0, := 0.0;
	var cordlga1 >= 0.0, := 0.0;
	var cordlga2 >= 0.0, := 0.0;
	var cordlga3 >= 0.0, := 0.0;
	var clgaord0 >= 0.0, := 0.0;
	var clgaord1 >= 0.0, := 0.0;
	var clgaord2 >= 0.0, := 0.0;
	var clgaord3 >= 0.0, := 0.0;
	var clgacle0 >= 0.0, := 0.0;
	var clgacle1 >= 0.0, := 0.0;
	var clgacle2 >= 0.0, := 0.0;
	var cclelga0 >= 0.0, := 0.0;
	var cclelga1 >= 0.0, := 0.0;
	var cclelga2 >= 0.0, := 0.0;
	var cclelga3 >= 0.0, := 0.0;
	var grdtimo1 >= 0.0, := 0.0;
	var grdtimn1 >= -100.0 ,  <= 0.0, := 0.0;
	var grdtimo2 >= 0.0, := 0.0;
	var grdtimn2 >= -90.0 ,  <= 0.0, := 0.0;
	var grdtimn3 >= -45.0 ,  <= 0.0, := 0.0;
	var grdtimn4 >= -45.0 ,  <= 0.0, := 0.0;
	var n1003ac1 >= 0.0, := 0.0 ,  <= 7.0;
	var n1003ac2 >= 0.0, := 0.0 ,  <= 7.0;
	var n1003ac3 >= 0.0, := 0.0 ,  <= 2.0;
	var n1003ac4 >= 0.0, := 0.0 ,  <= 2.0;
	var n1005ac1 >= 0.0, := 0.0 ,  <= 7.0;
	var n1005ac2 >= 0.0, := 0.0 ,  <= 7.0;
	var n1005ac3 >= 0.0, := 0.0 ,  <= 2.0;
	var n1005ac4 >= 0.0, := 0.0 ,  <= 2.0;
	var n1011ac1 >= 0.0, := 0.0 ,  <= 7.0;
	var n1011ac2 >= 0.0, := 0.0 ,  <= 7.0;
	var n1011ac3 >= 0.0, := 0.0 ,  <= 2.0;
	var n1011ac4 >= 0.0, := 0.0 ,  <= 2.0;
	var n1013ac2 >= 0.0, := 0.0 ,  <= 14.0;
	var n1013ac4 >= 0.0, := 0.0 ,  <= 2.0;
	var n1015ac2 >= 0.0, := 0.0 ,  <= 7.0;
	var n1015ac4 >= 0.0, := 0.0 ,  <= 2.0;
	var n1017ac2 >= 0.0, := 0.0 ,  <= 7.0;
	var n1017ac4 >= 0.0, := 0.0 ,  <= 2.0;
	var n1019ac2 >= 0.0, := 0.0 ,  <= 7.0;
	var n1019ac4 >= 0.0, := 0.0 ,  <= 2.0;
	var n1021ac1 >= 0.0, := 0.0 ,  <= 7.0;
	var n1021ac2 >= 0.0, := 0.0 ,  <= 7.0;
	var n1021ac3 >= 0.0, := 0.0 ,  <= 2.0;
	var n1021ac4 >= 0.0, := 0.0 ,  <= 2.0;
	var n1002ac1 >= 0.0, := 0.0 ,  <= 7.0;
	var n1002ac2 >= 0.0, := 0.0 ,  <= 7.0;
	var n1002ac3 >= 0.0, := 0.0 ,  <= 2.0;
	var n1002ac4 >= 0.0, := 0.0 ,  <= 2.0;
	var n1004ac2 >= 0.0, := 0.0 ,  <= 7.0;
	var n1004ac4 >= 0.0, := 0.0 ,  <= 2.0;
	var n1006ac1 >= 0.0, := 0.0 ,  <= 7.0;
	var n1006ac2 >= 0.0, := 0.0 ,  <= 7.0;
	var n1006ac3 >= 0.0, := 0.0 ,  <= 2.0;
	var n1006ac4 >= 0.0, := 0.0 ,  <= 2.0;
	var n1008ac2 >= 0.0, := 0.0 ,  <= 14.0;
	var n1008ac4 >= 0.0, := 0.0 ,  <= 2.0;
	var n1010ac2 >= 0.0, := 0.0 ,  <= 14.0;
	var n1010ac4 >= 0.0, := 0.0 ,  <= 2.0;
	var n1012ac2 >= 0.0, := 0.0 ,  <= 7.0;
	var n1012ac4 >= 0.0, := 0.0 ,  <= 2.0;
	var n1014ac2 >= 0.0, := 0.0 ,  <= 14.0;
	var n1014ac4 >= 0.0, := 0.0 ,  <= 2.0;
	var n1100ac2 >= 0.0, := 0.0 ,  <= 7.0;
	var n1100ac4 >= 0.0, := 0.0 ,  <= 7.0;
	var n1102ac2 >= 0.0, := 0.0 ,  <= 7.0;
	var n1102ac4 >= 0.0, := 0.0 ,  <= 7.0;
	var n1200ac2 >= 0.0, := 0.0 ,  <= 14.0;
	var n1200ac4 >= 0.0, := 0.0 ,  <= 7.0;
	var n1201ac2 >= 0.0, := 0.0 ,  <= 14.0;
	var n1201ac4 >= 0.0, := 0.0 ,  <= 7.0;

minimize obj:
	 - 1.0 * pbosord0 * pbosord0 + -0.922535241 * pbosord1 * pbosord1 + 
	-0.845070422 * pbosord2 * pbosord2 + -0.767605662 * pbosord3 * pbosord3 + 
	-0.690140843 * pbosord4 * pbosord4 + -0.612676084 * pboslga0 * pboslga0 + 
	-0.535211265 * pboslga1 * pboslga1 + -0.457746476 * pboslga2 * pboslga2 + 
	-0.380281687 * pboslga3 * pboslga3 + -0.302816898 * pboscle0 * pboscle0 + 
	-0.225352108 * pboscle1 * pboscle1 + -0.147887319 * pboscle2 * pboscle2 + 
	-0.070422538 * pboscle3 * pboscle3 + 0.007042253 * pordbos0 * pordbos0 + 
	0.084507041 * pordbos1 * pordbos1 + 0.161971837 * pordbos2 * pordbos2 + 
	0.239436626 * pordbos3 * pordbos3 + 0.316901416 * pordbos4 * pordbos4 + 
	0.394366205 * pordlga0 * pordlga0 + 0.471830994 * pordlga1 * pordlga1 + 
	0.549295783 * pordlga2 * pordlga2 + 0.626760542 * pordlga3 * pordlga3 + 
	0.704225361 * pordcle0 * pordcle0 + 0.781690121 * pordcle1 * pordcle1 + 
	0.85915494 * pordcle2 * pordcle2 + 0.936619699 * pordcle3 * pordcle3 + 
	1.014084458 * plgabos0 * plgabos0 + 1.091549277 * plgabos1 * plgabos1 + 
	1.169014096 * plgabos2 * plgabos2 + 1.246478915 * plgabos3 * plgabos3 + 
	1.323943615 * plgabos4 * plgabos4 + 1.401408434 * plgabos5 * plgabos5 + 
	1.478873253 * plgaord0 * plgaord0 + 1.556338072 * plgaord1 * plgaord1 + 
	1.633802772 * plgaord2 * plgaord2 + 1.711267591 * plgaord3 * plgaord3 + 
	1.788732409 * plgacle0 * plgacle0 + 1.866197228 * plgacle1 * plgacle1 + 
	1.943661928 * plgacle2 * plgacle2 + 2.021126747 * pclebos0 * pclebos0 + 
	2.098591566 * pclebos1 * pclebos1 + 2.176056385 * pclebos2 * pclebos2 + 
	2.253521204 * pclebos3 * pclebos3 + 2.330986023 * pclebos4 * pclebos4 + 
	2.408450603 * pcleord0 * pcleord0 + 2.485915422 * pcleord1 * pcleord1 + 
	2.563380241 * pcleord2 * pcleord2 + 2.64084506 * pcleord3 * pcleord3 + 
	2.718309879 * pcleord4 * pcleord4 + 2.795774698 * pclelga0 * pclelga0 + 
	2.873239517 * pclelga1 * pclelga1 + 2.950704336 * pclelga2 * pclelga2 + 
	3.028168917 * pclelga3 * pclelga3 + 3.105633736 * bbosord0 * bbosord0 + 
	3.183098555 * bbosord1 * bbosord1 + 3.260563374 * bboscle0 * bboscle0 + 
	3.338028193 * bboscle1 * bboscle1 + 3.415493011 * bordbos0 * bordbos0 + 
	3.49295783 * cbosord0 * cbosord0 + 3.570422649 * cbosord1 * cbosord1 + 
	3.64788723 * cbosord2 * cbosord2 + 3.725352049 * cbosord3 * cbosord3 + 
	3.802816868 * cbosord4 * cbosord4 + 3.880281687 * cboscle0 * cboscle0 + 
	3.957746506 * cboscle1 * cboscle1 + 4.035211086 * cboscle2 * cboscle2 + 
	4.112676144 * cboscle3 * cboscle3 + 4.190140724 * cordbos0 * cordbos0 + 
	4.267605782 * cordbos1 * cordbos1 + 4.345070362 * cordbos2 * cordbos2 + 
	4.422535419 * cordbos3 * cordbos3 + 4.5 * cordbos4 * cordbos4 + 4.577464581 * 
	cordlga0 * cordlga0 + 4.654929638 * cordlga1 * cordlga1 + 4.732394218 * 
	cordlga2 * cordlga2 + 4.809859276 * cordlga3 * cordlga3 + 4.887323856 * 
	clgaord0 * clgaord0 + 4.964788914 * clgaord1 * clgaord1 + 5.042253494 * 
	clgaord2 * clgaord2 + 5.119718075 * clgaord3 * clgaord3 + 5.197183132 * 
	clgacle0 * clgacle0 + 5.274647713 * clgacle1 * clgacle1 + 5.35211277 * clgacle2 
	* clgacle2 + 5.429577351 * cclelga0 * cclelga0 + 5.507042408 * cclelga1 * 
	cclelga1 + 5.584506989 * cclelga2 * cclelga2 + 5.661972046 * cclelga3 * 
	cclelga3 + 5.739436626 * grdtimo1 * grdtimo1 + 5.816901207 * grdtimn1 * 
	grdtimn1 + 5.894366264 * grdtimo2 * grdtimo2 + 5.971830845 * grdtimn2 * 
	grdtimn2 + 6.049295902 * grdtimn3 * grdtimn3 + 6.126760483 * grdtimn4 * 
	grdtimn4 + 6.20422554 * n1003ac1 * n1003ac1 + 6.281690121 * n1003ac2 * n1003ac2 
	+ 6.359154701 * n1003ac3 * n1003ac3 + 6.436619759 * n1003ac4 * n1003ac4 + 
	6.514084339 * n1005ac1 * n1005ac1 + 6.591549397 * n1005ac2 * n1005ac2 + 
	6.669013977 * n1005ac3 * n1005ac3 + 6.746479034 * n1005ac4 * n1005ac4 + 
	6.823943615 * n1011ac1 * n1011ac1 + 6.901408672 * n1011ac2 * n1011ac2 + 
	6.978873253 * n1011ac3 * n1011ac3 + 7.056337833 * n1011ac4 * n1011ac4 + 
	7.133802891 * n1013ac2 * n1013ac2 + 7.211267471 * n1013ac4 * n1013ac4 + 
	7.288732529 * n1015ac2 * n1015ac2 + 7.366197109 * n1015ac4 * n1015ac4 + 
	7.443662167 * n1017ac2 * n1017ac2 + 7.521126747 * n1017ac4 * n1017ac4 + 
	7.598591328 * n1019ac2 * n1019ac2 + 7.676056385 * n1019ac4 * n1019ac4 + 
	7.753520966 * n1021ac1 * n1021ac1 + 7.830986023 * n1021ac2 * n1021ac2 + 
	7.908450603 * n1021ac3 * n1021ac3 + 7.985915661 * n1021ac4 * n1021ac4 + 
	8.063380241 * n1002ac1 * n1002ac1 + 8.140845299 * n1002ac2 * n1002ac2 + 
	8.218309402 * n1002ac3 * n1002ac3 + 8.29577446 * n1002ac4 * n1002ac4 + 
	8.373239517 * n1004ac2 * n1004ac2 + 8.450704575 * n1004ac4 * n1004ac4 + 
	8.528168678 * n1006ac1 * n1006ac1 + 8.605633736 * n1006ac2 * n1006ac2 + 
	8.683098793 * n1006ac3 * n1006ac3 + 8.76056385 * n1006ac4 * n1006ac4 + 
	8.838027954 * n1008ac2 * n1008ac2 + 8.915493011 * n1008ac4 * n1008ac4 + 
	8.992958069 * n1010ac2 * n1010ac2 + 9.070422173 * n1010ac4 * n1010ac4 + 
	9.14788723 * n1012ac2 * n1012ac2 + 9.225352287 * n1012ac4 * n1012ac4 + 
	9.302817345 * n1014ac2 * n1014ac2 + 9.380281448 * n1014ac4 * n1014ac4 + 
	9.457746506 * n1100ac2 * n1100ac2 + 9.535211563 * n1100ac4 * n1100ac4 + 
	9.612675667 * n1102ac2 * n1102ac2 + 9.690140724 * n1102ac4 * n1102ac4 + 
	9.767605782 * n1200ac2 * n1200ac2 + 9.845070839 * n1200ac4 * n1200ac4 + 
	9.922534943 * n1201ac2 * n1201ac2 + 10.0 * n1201ac4 * n1201ac4 - 0.075*pbosord0 
	- 0.075*pbosord1 - 0.075*pbosord2 - 0.075*pbosord3 - 0.075*pbosord4 - 
	0.027*pboslga0 - 0.027*pboslga1 - 0.027*pboslga2 - 0.027*pboslga3 - 
	0.053*pboscle0 - 0.053*pboscle1 - 0.053*pboscle2 - 0.053*pboscle3 - 
	0.075*pordbos0 - 0.075*pordbos1 - 0.075*pordbos2 - 0.075*pordbos3 - 
	0.075*pordbos4 - 0.068*pordlga0 - 0.068*pordlga1 - 0.068*pordlga2 - 
	0.068*pordlga3 - 0.035*pordcle0 - 0.035*pordcle1 - 0.035*pordcle2 - 
	0.035*pordcle3 - 0.027*plgabos0 - 0.027*plgabos1 - 0.027*plgabos2 - 
	0.027*plgabos3 - 0.027*plgabos4 - 0.027*plgabos5 - 0.068*plgaord0 - 
	0.068*plgaord1 - 0.068*plgaord2 - 0.068*plgaord3 - 0.037*plgacle0 - 
	0.037*plgacle1 - 0.037*plgacle2 - 0.053*pclebos0 - 0.053*pclebos1 - 
	0.053*pclebos2 - 0.053*pclebos3 - 0.053*pclebos4 - 0.035*pcleord0 - 
	0.035*pcleord1 - 0.035*pcleord2 - 0.035*pcleord3 - 0.035*pcleord4 - 
	0.037*pclelga0 - 0.037*pclelga1 - 0.037*pclelga2 - 0.037*pclelga3 - 
	0.75*bbosord0 - 0.75*bbosord1 - 0.53*bboscle0 - 0.53*bboscle1 - 0.75*bordbos0 - 
	0.75*cbosord0 - 0.75*cbosord1 - 0.75*cbosord2 - 0.75*cbosord3 - 0.75*cbosord4 - 
	0.53*cboscle0 - 0.53*cboscle1 - 0.53*cboscle2 - 0.53*cboscle3 - 0.75*cordbos0 - 
	0.75*cordbos1 - 0.75*cordbos2 - 0.75*cordbos3 - 0.75*cordbos4 - 0.68*cordlga0 - 
	0.68*cordlga1 - 0.68*cordlga2 - 0.68*cordlga3 - 0.68*clgaord0 - 0.68*clgaord1 - 
	0.68*clgaord2 - 0.68*clgaord3 - 0.37*clgacle0 - 0.37*clgacle1 - 0.37*clgacle2 - 
	0.37*cclelga0 - 0.37*cclelga1 - 0.37*cclelga2 - 0.37*cclelga3 + 0.65*grdtimo1 - 
	0.05*grdtimn1 + 0.275*grdtimo2 - 0.02889*grdtimn2 - 0.03611*grdtimn3 - 
	0.01333*grdtimn4 + 0.01*n1003ac1 + 0.02*n1003ac2 + 0.03*n1003ac3 + 
	0.04*n1003ac4 + 4.04337*n1005ac1 + 1.60964*n1005ac2 + 1.93119*n1005ac3 + 
	1.34618*n1005ac4 + 1.1*n1011ac1 + 1.2*n1011ac2 + 1.3*n1011ac3 + 
	2.34647*n1011ac4 + 2.36783*n1013ac2 + 1.90292*n1013ac4 + 1.36416*n1015ac2 + 
	1.14401*n1015ac4 + 1.64736*n1017ac2 + 1.33312*n1017ac4 + 2.55381*n1019ac2 + 
	2.14431*n1019ac4 + 3.12679*n1021ac1 + 1.18965*n1021ac2 + 1.42472*n1021ac3 + 
	1.0003*n1021ac4 + 3.12679*n1002ac1 + 1.18965*n1002ac2 + 1.42472*n1002ac3 + 
	1.0003*n1002ac4 + 1.64736*n1004ac2 + 1.33312*n1004ac4 + 7.17016*n1006ac1 + 
	2.79929*n1006ac2 + 3.35591*n1006ac3 + 2.34647*n1006ac4 + 2.36783*n1008ac2 + 
	1.90292*n1008ac4 + 3.27428*n1010ac2 + 2.71411*n1010ac4 + 1.36416*n1012ac2 + 
	1.14401*n1012ac4 + 2.08463*n1014ac2 + 1.71382*n1014ac4 + 6.07357*n1100ac2 + 
	5.06059*n1100ac4 + 5.16712*n1102ac2 + 4.24939*n1102ac4 + 0.72047*n1200ac2 + 
	0.5698*n1200ac4 + 0.72047*n1201ac2 + 0.5698*n1201ac4;

subject to revenues:
	0 <= 0.075*pbosord0 + 0.075*pbosord1 + 0.075*pbosord2 + 0.075*pbosord3 + 
	0.075*pbosord4 + 0.027*pboslga0 + 0.027*pboslga1 + 0.027*pboslga2 + 
	0.027*pboslga3 + 0.053*pboscle0 + 0.053*pboscle1 + 0.053*pboscle2 + 
	0.053*pboscle3 + 0.075*pordbos0 + 0.075*pordbos1 + 0.075*pordbos2 + 
	0.075*pordbos3 + 0.075*pordbos4 + 0.068*pordlga0 + 0.068*pordlga1 + 
	0.068*pordlga2 + 0.068*pordlga3 + 0.035*pordcle0 + 0.035*pordcle1 + 
	0.035*pordcle2 + 0.035*pordcle3 + 0.027*plgabos0 + 0.027*plgabos1 + 
	0.027*plgabos2 + 0.027*plgabos3 + 0.027*plgabos4 + 0.027*plgabos5 + 
	0.068*plgaord0 + 0.068*plgaord1 + 0.068*plgaord2 + 0.068*plgaord3 + 
	0.037*plgacle0 + 0.037*plgacle1 + 0.037*plgacle2 + 0.053*pclebos0 + 
	0.053*pclebos1 + 0.053*pclebos2 + 0.053*pclebos3 + 0.053*pclebos4 + 
	0.035*pcleord0 + 0.035*pcleord1 + 0.035*pcleord2 + 0.035*pcleord3 + 
	0.035*pcleord4 + 0.037*pclelga0 + 0.037*pclelga1 + 0.037*pclelga2 + 
	0.037*pclelga3 + 0.75*bbosord0 + 0.75*bbosord1 + 0.53*bboscle0 + 0.53*bboscle1 
	+ 0.75*bordbos0 + 0.75*cbosord0 + 0.75*cbosord1 + 0.75*cbosord2 + 0.75*cbosord3 
	+ 0.75*cbosord4 + 0.53*cboscle0 + 0.53*cboscle1 + 0.53*cboscle2 + 0.53*cboscle3 
	+ 0.75*cordbos0 + 0.75*cordbos1 + 0.75*cordbos2 + 0.75*cordbos3 + 0.75*cordbos4 
	+ 0.68*cordlga0 + 0.68*cordlga1 + 0.68*cordlga2 + 0.68*cordlga3 + 0.68*clgaord0 
	+ 0.68*clgaord1 + 0.68*clgaord2 + 0.68*clgaord3 + 0.37*clgacle0 + 0.37*clgacle1 
	+ 0.37*clgacle2 + 0.37*cclelga0 + 0.37*cclelga1 + 0.37*cclelga2 + 0.37*cclelga3;
subject to acocosts:
	0 <= 0.65*grdtimo1 - 0.05*grdtimn1 + 0.275*grdtimo2 - 0.02889*grdtimn2 - 
	0.03611*grdtimn3 - 0.01333*grdtimn4 + 0.01*n1003ac1 + 0.02*n1003ac2 + 
	0.03*n1003ac3 + 0.04*n1003ac4 + 4.04337*n1005ac1 + 1.60964*n1005ac2 + 
	1.93119*n1005ac3 + 1.34618*n1005ac4 + 1.1*n1011ac1 + 1.2*n1011ac2 + 
	1.3*n1011ac3 + 2.34647*n1011ac4 + 2.36783*n1013ac2 + 1.90292*n1013ac4 + 
	1.36416*n1015ac2 + 1.14401*n1015ac4 + 1.64736*n1017ac2 + 1.33312*n1017ac4 + 
	2.55381*n1019ac2 + 2.14431*n1019ac4 + 3.12679*n1021ac1 + 1.18965*n1021ac2 + 
	1.42472*n1021ac3 + 1.0003*n1021ac4 + 3.12679*n1002ac1 + 1.18965*n1002ac2 + 
	1.42472*n1002ac3 + 1.0003*n1002ac4 + 1.64736*n1004ac2 + 1.33312*n1004ac4 + 
	7.17016*n1006ac1 + 2.79929*n1006ac2 + 3.35591*n1006ac3 + 2.34647*n1006ac4 + 
	2.36783*n1008ac2 + 1.90292*n1008ac4 + 3.27428*n1010ac2 + 2.71411*n1010ac4 + 
	1.36416*n1012ac2 + 1.14401*n1012ac4 + 2.08463*n1014ac2 + 1.71382*n1014ac4 + 
	6.07357*n1100ac2 + 5.06059*n1100ac4 + 5.16712*n1102ac2 + 4.24939*n1102ac4 + 
	0.72047*n1200ac2 + 0.5698*n1200ac4 + 0.72047*n1201ac2 + 0.5698*n1201ac4;
subject to fuelaval:
	0 >= 7.98429*n1003ac1 + 2.51914*n1003ac2 + 3.15178*n1003ac3 + 1.99337*n1003ac4 
	+ 5.83404*n1005ac1 + 1.82258*n1005ac2 + 2.27351*n1005ac3 + 1.41795*n1005ac4 + 
	9.91398*n1011ac1 + 3.07692*n1011ac2 + 3.83055*n1011ac3 + 2.3665*n1011ac4 + 
	3.16965*n1013ac2 + 2.4431*n1013ac4 + 1.49045*n1015ac2 + 1.14359*n1015ac4 + 
	2.21183*n1017ac2 + 1.73951*n1017ac4 + 2.7448*n1019ac2 + 2.09214*n1019ac4 + 
	4.07994*n1021ac1 + 1.25435*n1021ac2 + 1.55704*n1021ac3 + 0.94855*n1021ac4 + 
	4.07994*n1002ac1 + 1.25435*n1002ac2 + 1.55704*n1002ac3 + 0.94855*n1002ac4 + 
	2.21183*n1004ac2 + 1.73951*n1004ac4 + 9.91398*n1006ac1 + 3.07692*n1006ac2 + 
	3.83055*n1006ac3 + 2.3665*n1006ac4 + 3.16965*n1008ac2 + 2.4431*n1008ac4 + 
	3.70262*n1010ac2 + 2.79573*n1010ac4 + 1.49045*n1012ac2 + 1.14359*n1012ac4 + 
	2.44827*n1014ac2 + 1.84718*n1014ac4 + 6.77953*n1100ac2 + 5.16223*n1100ac4 + 
	6.24657*n1102ac2 + 4.8096*n1102ac4 + 0.95782*n1200ac2 + 0.70359*n1200ac4 + 
	0.95782*n1201ac2 + 0.70359*n1201ac4 + 100000.0*-1;
subject to systdept:
	0 <= n1003ac1 + n1003ac2 + n1003ac3 + n1003ac4 + n1005ac1 + n1005ac2 + n1005ac3 
	+ n1005ac4 + 2.0*n1011ac1 + 2.0*n1011ac2 + 2.0*n1011ac3 + 2.0*n1011ac4 + 
	2.0*n1013ac2 + 2.0*n1013ac4 + n1015ac2 + n1015ac4 + n1017ac2 + n1017ac4 + 
	2.0*n1019ac2 + 2.0*n1019ac4 + n1021ac1 + n1021ac2 + n1021ac3 + n1021ac4 + 
	n1002ac1 + n1002ac2 + n1002ac3 + n1002ac4 + n1004ac2 + n1004ac4 + 2.0*n1006ac1 
	+ 2.0*n1006ac2 + 2.0*n1006ac3 + 2.0*n1006ac4 + 2.0*n1008ac2 + 2.0*n1008ac4 + 
	3.0*n1010ac2 + 3.0*n1010ac4 + n1012ac2 + n1012ac4 + 2.0*n1014ac2 + 2.0*n1014ac4 
	+ 5.0*n1100ac2 + 5.0*n1100ac4 + 4.0*n1102ac2 + 4.0*n1102ac4 + n1200ac2 + 
	n1200ac4 + n1201ac2 + n1201ac4 + 50.0*-1;
subject to acmiles:
	0 <= 0.86441*n1003ac1 + 0.86441*n1003ac2 + 0.86441*n1003ac3 + 0.86441*n1003ac4 
	+ 0.56156*n1005ac1 + 0.56156*n1005ac2 + 0.56156*n1005ac3 + 0.56156*n1005ac4 + 
	0.87605*n1011ac1 + 0.87605*n1011ac2 + 0.87605*n1011ac3 + 0.87605*n1011ac4 + 
	0.91637*n1013ac2 + 0.91637*n1013ac4 + 0.41715*n1015ac2 + 0.41715*n1015ac4 + 
	0.7308*n1017ac2 + 0.7308*n1017ac4 + 0.73165*n1019ac2 + 0.73165*n1019ac4 + 
	0.3145*n1021ac1 + 0.3145*n1021ac2 + 0.3145*n1021ac3 + 0.3145*n1021ac4 + 
	0.3145*n1002ac1 + 0.3145*n1002ac2 + 0.3145*n1002ac3 + 0.3145*n1002ac4 + 
	0.7308*n1004ac2 + 0.7308*n1004ac4 + 0.87605*n1006ac1 + 0.87605*n1006ac2 + 
	0.87605*n1006ac3 + 0.87605*n1006ac4 + 0.91637*n1008ac2 + 0.91637*n1008ac4 + 
	0.91722*n1010ac2 + 0.91722*n1010ac4 + 0.41715*n1012ac2 + 0.41715*n1012ac4 + 
	0.60273*n1014ac2 + 0.60273*n1014ac4 + 1.79328*n1100ac2 + 1.79328*n1100ac4 + 
	1.79242*n1102ac2 + 1.79242*n1102ac4 + 0.18557*n1200ac2 + 0.18557*n1200ac4 + 
	0.18557*n1201ac2 + 0.18557*n1201ac4;
subject to asmiles:
	0 <= 351.81396*n1003ac1 + 113.23743*n1003ac2 + 228.55299*n1005ac1 + 
	73.56374*n1005ac2 + 356.55371*n1011ac1 + 114.76299*n1011ac2 + 
	120.04449*n1013ac2 + 54.64705*n1015ac2 + 95.73444*n1017ac2 + 95.8463*n1019ac2 + 
	128.00075*n1021ac1 + 41.19926*n1021ac2 + 128.00075*n1002ac1 + 41.19926*n1002ac2 
	+ 95.73444*n1004ac2 + 356.55371*n1006ac1 + 114.76299*n1006ac2 + 
	120.04449*n1008ac2 + 120.15637*n1010ac2 + 54.64705*n1012ac2 + 78.95706*n1014ac2 
	+ 234.91937*n1100ac2 + 234.80756*n1102ac2 + 24.31007*n1200ac2 + 
	24.31007*n1201ac2;
subject to passngrs:
	0 <= pbosord0 + pbosord1 + pbosord2 + pbosord3 + pbosord4 + pboslga0 + pboslga1 
	+ pboslga2 + pboslga3 + pboscle0 + pboscle1 + pboscle2 + pboscle3 + pordbos0 + 
	pordbos1 + pordbos2 + pordbos3 + pordbos4 + pordlga0 + pordlga1 + pordlga2 + 
	pordlga3 + pordcle0 + pordcle1 + pordcle2 + pordcle3 + plgabos0 + plgabos1 + 
	plgabos2 + plgabos3 + plgabos4 + plgabos5 + plgaord0 + plgaord1 + plgaord2 + 
	plgaord3 + plgacle0 + plgacle1 + plgacle2 + pclebos0 + pclebos1 + pclebos2 + 
	pclebos3 + pclebos4 + pcleord0 + pcleord1 + pcleord2 + pcleord3 + pcleord4 + 
	pclelga0 + pclelga1 + pclelga2 + pclelga3 + 9431.0*-1;
subject to rpmiles:
	0 <= 0.86441*pbosord0 + 0.87605*pbosord1 + 0.91637*pbosord2 + 0.91722*pbosord3 
	+ 0.87605*pbosord4 + 0.18557*pboslga0 + 0.18557*pboslga1 + 1.60685*pboslga2 + 
	0.18557*pboslga3 + 0.56156*pboscle0 + 0.56156*pboscle1 + 0.60273*pboscle2 + 
	0.56156*pboscle3 + 0.87605*pordbos0 + 0.91637*pordbos1 + 0.91722*pordbos2 + 
	0.87605*pordbos3 + 0.91637*pordbos4 + 0.7308*pordlga0 + 0.7308*pordlga1 + 
	0.73165*pordlga2 + 0.7308*pordlga3 + 0.3145*pordcle0 + 0.3145*pordcle1 + 
	0.3145*pordcle2 + 0.3145*pordcle3 + 0.18557*plgabos0 + 0.18557*plgabos1 + 
	0.18557*plgabos2 + 1.6077*plgabos3 + 0.18557*plgabos4 + 0.18557*plgabos5 + 
	0.7308*plgaord0 + 0.7308*plgaord1 + 0.73165*plgaord2 + 0.73165*plgaord3 + 
	0.41715*plgacle0 + 0.41715*plgacle1 + 0.41715*plgacle2 + 0.56156*pclebos0 + 
	0.60273*pclebos1 + 0.60273*pclebos2 + 0.56156*pclebos3 + 1.23087*pclebos4 + 
	0.3145*pcleord0 + 0.3145*pcleord1 + 0.3145*pcleord2 + 0.3145*pcleord3 + 
	0.3145*pcleord4 + 0.41715*pclelga0 + 0.41715*pclelga1 + 0.41715*pclelga2 + 
	1.0453*pclelga3;
subject to lfrpmasm:
	0 <= -0.86441*pbosord0 - 0.87605*pbosord1 - 0.91637*pbosord2 - 0.91722*pbosord3 
	- 0.87605*pbosord4 - 0.18557*pboslga0 - 0.18557*pboslga1 - 1.60685*pboslga2 - 
	0.18557*pboslga3 - 0.56156*pboscle0 - 0.56156*pboscle1 - 0.60273*pboscle2 - 
	0.56156*pboscle3 - 0.87605*pordbos0 - 0.91637*pordbos1 - 0.91722*pordbos2 - 
	0.87605*pordbos3 - 0.91637*pordbos4 - 0.7308*pordlga0 - 0.7308*pordlga1 - 
	0.73165*pordlga2 - 0.7308*pordlga3 - 0.3145*pordcle0 - 0.3145*pordcle1 - 
	0.3145*pordcle2 - 0.3145*pordcle3 - 0.18557*plgabos0 - 0.18557*plgabos1 - 
	0.18557*plgabos2 - 1.6077*plgabos3 - 0.18557*plgabos4 - 0.18557*plgabos5 - 
	0.7308*plgaord0 - 0.7308*plgaord1 - 0.73165*plgaord2 - 0.73165*plgaord3 - 
	0.41715*plgacle0 - 0.41715*plgacle1 - 0.41715*plgacle2 - 0.56156*pclebos0 - 
	0.60273*pclebos1 - 0.60273*pclebos2 - 0.56156*pclebos3 - 1.23087*pclebos4 - 
	0.3145*pcleord0 - 0.3145*pcleord1 - 0.3145*pcleord2 - 0.3145*pcleord3 - 
	0.3145*pcleord4 - 0.41715*pclelga0 - 0.41715*pclelga1 - 0.41715*pclelga2 - 
	1.0453*pclelga3 + 211.088376*n1003ac1 + 67.942458*n1003ac2 + 
	137.131794*n1005ac1 + 44.138244*n1005ac2 + 213.93222599999999*n1011ac1 + 
	68.857794*n1011ac2 + 72.026694*n1013ac2 + 32.78823*n1015ac2 + 
	57.440664*n1017ac2 + 57.50778*n1019ac2 + 76.80045*n1021ac1 + 24.719556*n1021ac2 
	+ 76.80045*n1002ac1 + 24.719556*n1002ac2 + 57.440664*n1004ac2 + 
	213.93222599999999*n1006ac1 + 68.857794*n1006ac2 + 72.026694*n1008ac2 + 
	72.093822*n1010ac2 + 32.78823*n1012ac2 + 47.374236*n1014ac2 + 
	140.951622*n1100ac2 + 140.884536*n1102ac2 + 14.586042*n1200ac2 + 
	14.586042*n1201ac2;
subject to atonmile:
	0 <= 25.93224*n1003ac1 + 34.57631*n1003ac3 + 12.96612*n1003ac4 + 
	16.84665*n1005ac1 + 22.4622*n1005ac3 + 8.42333*n1005ac4 + 26.2816*n1011ac1 + 
	35.04214*n1011ac3 + 13.1408*n1011ac4 + 13.74556*n1013ac4 + 6.25729*n1015ac4 + 
	10.96196*n1017ac4 + 10.97477*n1019ac4 + 9.43495*n1021ac1 + 12.57993*n1021ac3 + 
	4.71747*n1021ac4 + 9.43495*n1002ac1 + 12.57993*n1002ac3 + 4.71747*n1002ac4 + 
	10.96196*n1004ac4 + 26.2816*n1006ac1 + 35.04214*n1006ac3 + 13.1408*n1006ac4 + 
	13.74556*n1008ac4 + 13.75836*n1010ac4 + 6.25729*n1012ac4 + 9.04089*n1014ac4 + 
	26.89915*n1100ac4 + 26.88635*n1102ac4 + 2.7836*n1200ac4 + 2.7836*n1201ac4;
subject to rtonmile:
	0 <= 0.86441*bbosord0 + 0.87605*bbosord1 + 0.56156*bboscle0 + 0.56156*bboscle1 
	+ 0.87605*bordbos0 + 0.86441*cbosord0 + 0.87605*cbosord1 + 0.91637*cbosord2 + 
	0.91722*cbosord3 + 0.87605*cbosord4 + 0.56156*cboscle0 + 0.56156*cboscle1 + 
	0.60273*cboscle2 + 0.56156*cboscle3 + 0.87605*cordbos0 + 0.91637*cordbos1 + 
	0.91722*cordbos2 + 0.87605*cordbos3 + 0.91637*cordbos4 + 0.7308*cordlga0 + 
	0.7308*cordlga1 + 0.73165*cordlga2 + 0.7308*cordlga3 + 0.7308*clgaord0 + 
	0.7308*clgaord1 + 0.73165*clgaord2 + 0.73165*clgaord3 + 0.41715*clgacle0 + 
	0.41715*clgacle1 + 0.41715*clgacle2 + 0.41715*cclelga0 + 0.41715*cclelga1 + 
	0.41715*cclelga2 + 1.0453*cclelga3;
subject to lftnmile:
	0 <= -0.86441*bbosord0 - 0.87605*bbosord1 - 0.56156*bboscle0 - 0.56156*bboscle1 
	- 0.87605*bordbos0 - 0.86441*cbosord0 - 0.87605*cbosord1 - 0.91637*cbosord2 - 
	0.91722*cbosord3 - 0.87605*cbosord4 - 0.56156*cboscle0 - 0.56156*cboscle1 - 
	0.60273*cboscle2 - 0.56156*cboscle3 - 0.87605*cordbos0 - 0.91637*cordbos1 - 
	0.91722*cordbos2 - 0.87605*cordbos3 - 0.91637*cordbos4 - 0.7308*cordlga0 - 
	0.7308*cordlga1 - 0.73165*cordlga2 - 0.7308*cordlga3 - 0.7308*clgaord0 - 
	0.7308*clgaord1 - 0.73165*clgaord2 - 0.73165*clgaord3 - 0.41715*clgacle0 - 
	0.41715*clgacle1 - 0.41715*clgacle2 - 0.41715*cclelga0 - 0.41715*cclelga1 - 
	0.41715*cclelga2 - 1.0453*cclelga3 + 12.96612*n1003ac1 + 17.288155*n1003ac3 + 
	6.48306*n1003ac4 + 8.423325*n1005ac1 + 11.2311*n1005ac3 + 4.211665*n1005ac4 + 
	13.1408*n1011ac1 + 17.52107*n1011ac3 + 6.5704*n1011ac4 + 6.87278*n1013ac4 + 
	3.128645*n1015ac4 + 5.48098*n1017ac4 + 5.487385*n1019ac4 + 4.717475*n1021ac1 + 
	6.289965*n1021ac3 + 2.358735*n1021ac4 + 4.717475*n1002ac1 + 6.289965*n1002ac3 + 
	2.358735*n1002ac4 + 5.48098*n1004ac4 + 13.1408*n1006ac1 + 17.52107*n1006ac3 + 
	6.5704*n1006ac4 + 6.87278*n1008ac4 + 6.87918*n1010ac4 + 3.128645*n1012ac4 + 
	4.520445*n1014ac4 + 13.449575*n1100ac4 + 13.443175*n1102ac4 + 1.3918*n1200ac4 + 
	1.3918*n1201ac4;
subject to flav1:
	grdtimo1 + grdtimn1 + 2.32729*n1003ac1 + 1000.0*n1005ac1 + 25.0*n1011ac1 + 
	1.34295*n1021ac1 + 1.34295*n1002ac1 + 3.12813*n1006ac1 + 30.0*-1 = 0;
subject to flav2:
	grdtimo2 + grdtimn2 + 2.10966*n1003ac2 + 2000.0*n1005ac2 + 25.0*n1011ac2 + 
	2.65943*n1013ac2 + 1.25093*n1015ac2 + 1.85313*n1017ac2 + 2.30477*n1019ac2 + 
	1.05384*n1021ac2 + 1.05384*n1002ac2 + 1.85313*n1004ac2 + 2.58202*n1006ac2 + 
	2.65943*n1008ac2 + 3.11107*n1010ac2 + 1.25093*n1012ac2 + 2.05723*n1014ac2 + 
	5.69309*n1100ac2 + 5.24145*n1102ac2 + 0.8063*n1200ac2 + 0.8063*n1201ac2 + 
	45.0*-1 = 0;
subject to flav3:
	grdtimn3 + 2.10966*n1003ac3 + 3000.0*n1005ac3 + 26.0*n1011ac3 + 
	1.05384*n1021ac3 + 1.05384*n1002ac3 + 2.58202*n1006ac3 = 0;
subject to flav4:
	grdtimn4 + 2.10966*n1003ac4 + 1500.0*n1005ac4 + 27.0*n1011ac4 + 
	2.65943*n1013ac4 + 1.25093*n1015ac4 + 1.85313*n1017ac4 + 2.30477*n1019ac4 + 
	1.05384*n1021ac4 + 1.05384*n1002ac4 + 1.85313*n1004ac4 + 2.58202*n1006ac4 + 
	2.65943*n1008ac4 + 3.11107*n1010ac4 + 1.25093*n1012ac4 + 2.05723*n1014ac4 + 
	5.69309*n1100ac4 + 5.24145*n1102ac4 + 0.8063*n1200ac4 + 0.8063*n1201ac4 = 0;
subject to lf1003s1:
	0 <= -pbosord0 + 305.0*n1003ac1 + 98.0*n1003ac2;
subject to lf1003b1:
	0 <= -bbosord0 + 12.0*n1003ac1;
subject to lf1003c1:
	0 <= -cbosord0 + 30.0*n1003ac3 + 11.25*n1003ac4;
subject to lf1005s1:
	0 <= -pboscle0 + 305.0*n1005ac1 + 98.0*n1005ac2;
subject to lf1005b1:
	0 <= -bboscle0 + 12.0*n1005ac1;
subject to lf1005c1:
	0 <= -cboscle0 + 30.0*n1005ac3 + 11.25*n1005ac4;
subject to lf1011s1:
	0 <= -pbosord1 - pboscle1 + 305.0*n1011ac1 + 98.0*n1011ac2;
subject to lf1011s2:
	0 <= -pbosord1 - pcleord0 + 305.0*n1011ac1 + 98.0*n1011ac2;
subject to lf1011b1:
	0 <= -bbosord1 - bboscle1 + 12.0*n1011ac1;
subject to lf1011b2:
	0 <= -bbosord1 + 12.0*n1011ac1;
subject to lf1011c1:
	0 <= -cbosord1 - cboscle1 + 30.0*n1011ac3 + 11.25*n1011ac4;
subject to lf1011c2:
	0 <= -cbosord1 + 30.0*n1011ac3 + 11.25*n1011ac4;
subject to lf1013s1:
	0 <= -pbosord2 - pboslga0 + 98.0*n1013ac2;
subject to lf1013s2:
	0 <= -pbosord2 - plgaord0 + 98.0*n1013ac2;
subject to lf1013c1:
	0 <= -cbosord2 + 11.25*n1013ac4;
subject to lf1013c2:
	0 <= -cbosord2 - clgaord0 + 11.25*n1013ac4;
subject to lf1015s1:
	0 <= -plgacle0 + 98.0*n1015ac2;
subject to lf1015c1:
	0 <= -clgacle0 + 11.25*n1015ac4;
subject to lf1017s1:
	0 <= -plgaord1 + 98.0*n1017ac2;
subject to lf1017c1:
	0 <= -clgaord1 + 11.25*n1017ac4;
subject to lf1019s1:
	0 <= -plgaord2 - plgacle1 + 98.0*n1019ac2;
subject to lf1019s2:
	0 <= -plgaord2 - pcleord1 + 98.0*n1019ac2;
subject to lf1019c1:
	0 <= -clgaord2 - clgacle1 + 11.25*n1019ac4;
subject to lf1019c2:
	0 <= -clgaord2 + 11.25*n1019ac4;
subject to lf1021s1:
	0 <= -pcleord2 + 305.0*n1021ac1 + 98.0*n1021ac2;
subject to lf1021b1:
	0 <= 12.0*n1021ac1;
subject to lf1021c1:
	0 <= 30.0*n1021ac3 + 11.25*n1021ac4;
subject to lf1002s1:
	0 <= -pordcle0 + 305.0*n1002ac1 + 98.0*n1002ac2;
subject to lf1002b1:
	0 <= 12.0*n1002ac1;
subject to lf1002c1:
	0 <= 30.0*n1002ac3 + 11.25*n1002ac4;
subject to lf1004s1:
	0 <= -pordlga0 + 98.0*n1004ac2;
subject to lf1004c1:
	0 <= -cordlga0 + 11.25*n1004ac4;
subject to lf1006s1:
	0 <= -pordbos0 - pordcle1 + 305.0*n1006ac1 + 98.0*n1006ac2;
subject to lf1006s2:
	0 <= -pordbos0 - pclebos0 + 305.0*n1006ac1 + 98.0*n1006ac2;
subject to lf1006b1:
	0 <= -bordbos0 + 12.0*n1006ac1;
subject to lf1006b2:
	0 <= -bordbos0 + 12.0*n1006ac1;
subject to lf1006c1:
	0 <= -cordbos0 + 30.0*n1006ac3 + 11.25*n1006ac4;
subject to lf1006c2:
	0 <= -cordbos0 + 30.0*n1006ac3 + 11.25*n1006ac4;
subject to lf1008s1:
	0 <= -pordbos1 - pordlga1 + 98.0*n1008ac2;
subject to lf1008s2:
	0 <= -pordbos1 - plgabos0 + 98.0*n1008ac2;
subject to lf1008c1:
	0 <= -cordbos1 - cordlga1 + 11.25*n1008ac4;
subject to lf1008c2:
	0 <= -cordbos1 + 11.25*n1008ac4;
subject to lf1010s1:
	0 <= -pordbos2 - pordlga2 - pordcle2 + 98.0*n1010ac2;
subject to lf1010s2:
	0 <= -pordbos2 - pordlga2 - pclebos1 - pclelga0 + 98.0*n1010ac2;
subject to lf1010s3:
	0 <= -pordbos2 - plgabos1 - pclebos1 + 98.0*n1010ac2;
subject to lf1010c1:
	0 <= -cordbos2 - cordlga2 + 11.25*n1010ac4;
subject to lf1010c2:
	0 <= -cordbos2 - cordlga2 - cclelga0 + 11.25*n1010ac4;
subject to lf1010c3:
	0 <= -cordbos2 + 11.25*n1010ac4;
subject to lf1012s1:
	0 <= -pclelga1 + 98.0*n1012ac2;
subject to lf1012c1:
	0 <= -cclelga1 + 11.25*n1012ac4;
subject to lf1014s1:
	0 <= -pclebos2 - pclelga2 + 98.0*n1014ac2;
subject to lf1014s2:
	0 <= -plgabos2 - pclebos2 + 98.0*n1014ac2;
subject to lf1014c1:
	0 <= -cclelga2 + 11.25*n1014ac4;
subject to lf1014c2:
	0 <= 11.25*n1014ac4;
subject to lf1100s1:
	0 <= -pbosord3 - pboslga1 - pboscle2 + 98.0*n1100ac2;
subject to lf1100s2:
	0 <= -pbosord3 - pboscle2 - plgabos3 - plgaord3 - plgacle2 + 98.0*n1100ac2;
subject to lf1100s3:
	0 <= -pbosord3 - plgabos3 - plgaord3 - pcleord3 + 98.0*n1100ac2;
subject to lf1100s4:
	0 <= -pordbos3 - pordcle3 - plgabos3 + 98.0*n1100ac2;
subject to lf1100s5:
	0 <= -pordbos3 - plgabos3 - pclebos3 + 98.0*n1100ac2;
subject to lf1100c1:
	0 <= -cbosord3 - cboscle2 + 11.25*n1100ac4;
subject to lf1100c2:
	0 <= -cbosord3 - cboscle2 - clgaord3 - clgacle2 + 11.25*n1100ac4;
subject to lf1100c3:
	0 <= -cbosord3 - clgaord3 + 11.25*n1100ac4;
subject to lf1100c4:
	0 <= -cordbos3 + 11.25*n1100ac4;
subject to lf1100c5:
	0 <= -cordbos3 + 11.25*n1100ac4;
subject to lf1102s1:
	0 <= -pbosord4 - pboslga2 - pboscle3 + 98.0*n1102ac2;
subject to lf1102s2:
	0 <= -pbosord4 - pboslga2 - pclebos4 - pcleord4 - pclelga3 + 98.0*n1102ac2;
subject to lf1102s3:
	0 <= -pboslga2 - pordbos4 - pordlga3 - pclebos4 - pclelga3 + 98.0*n1102ac2;
subject to lf1102s4:
	0 <= -pordbos4 - plgabos4 - pclebos4 + 98.0*n1102ac2;
subject to lf1102c1:
	0 <= -cbosord4 - cboscle3 + 11.25*n1102ac4;
subject to lf1102c2:
	0 <= -cbosord4 - cclelga3 + 11.25*n1102ac4;
subject to lf1102c3:
	0 <= -cordbos4 - cordlga3 - cclelga3 + 11.25*n1102ac4;
subject to lf1102c4:
	0 <= -cordbos4 + 11.25*n1102ac4;
subject to lf1200s1:
	0 <= -pboslga3 + 98.0*n1200ac2;
subject to lf1200c1:
	0 <= 11.25*n1200ac4;
subject to lf1201s1:
	0 <= -plgabos5 + 98.0*n1201ac2;
subject to lf1201c1:
	0 <= 11.25*n1201ac4;
subject to noptcle0:
	0 <= n1005ac1 + n1005ac2 + n1005ac3 + n1005ac4 + 2.0*n1011ac1 + 2.0*n1011ac2 + 
	2.0*n1011ac3 + 2.0*n1011ac4 + n1015ac2 + n1015ac4 + 2.0*n1019ac2 + 2.0*n1019ac4 
	+ n1021ac1 + n1021ac2 + n1021ac3 + n1021ac4 + n1002ac1 + n1002ac2 + n1002ac3 + 
	n1002ac4 + 2.0*n1006ac1 + 2.0*n1006ac2 + 2.0*n1006ac3 + 2.0*n1006ac4 + 
	2.0*n1010ac2 + 2.0*n1010ac4 + n1012ac2 + n1012ac4 + n1014ac2 + n1014ac4 + 
	4.0*n1100ac2 + 4.0*n1100ac4 + 2.0*n1102ac2 + 2.0*n1102ac4 + 24.0*-1;
subject to contbos1:
	0 <= n1003ac1 + n1005ac1 + n1011ac1 - n1006ac1;
subject to contbos2:
	0 <= n1003ac2 + n1005ac2 + n1011ac2 + n1013ac2 - n1006ac2 - n1008ac2 - n1010ac2 
	- n1014ac2 + n1200ac2 - n1201ac2;
subject to contbos3:
	0 <= n1003ac3 + n1005ac3 + n1011ac3 - n1006ac3;
subject to contbos4:
	0 <= n1003ac4 + n1005ac4 + n1011ac4 + n1013ac4 - n1006ac4 - n1008ac4 - n1010ac4 
	- n1014ac4 + n1200ac4 - n1201ac4;
subject to contord1:
	0 <= -n1003ac1 - n1011ac1 - n1021ac1 + n1002ac1 + n1006ac1;
subject to contord2:
	0 <= -n1003ac2 - n1011ac2 - n1013ac2 - n1017ac2 - n1019ac2 - n1021ac2 + 
	n1002ac2 + n1004ac2 + n1006ac2 + n1008ac2 + n1010ac2;
subject to contord3:
	0 <= -n1003ac3 - n1011ac3 - n1021ac3 + n1002ac3 + n1006ac3;
subject to contord4:
	0 <= -n1003ac4 - n1011ac4 - n1013ac4 - n1017ac4 - n1019ac4 - n1021ac4 + 
	n1002ac4 + n1004ac4 + n1006ac4 + n1008ac4 + n1010ac4;
subject to contlga2:
	0 <= n1015ac2 + n1017ac2 + n1019ac2 - n1004ac2 - n1012ac2 - n1200ac2 + n1201ac2;
subject to contlga4:
	0 <= n1015ac4 + n1017ac4 + n1019ac4 - n1004ac4 - n1012ac4 - n1200ac4 + n1201ac4;
subject to contcle1:
	0 <= -n1005ac1 + n1021ac1 - n1002ac1;
subject to contcle2:
	0 <= -n1005ac2 - n1015ac2 + n1021ac2 - n1002ac2 + n1012ac2 + n1014ac2;
subject to contcle3:
	0 <= -n1005ac3 + n1021ac3 - n1002ac3;
subject to contcle4:
	0 <= -n1005ac4 - n1015ac4 + n1021ac4 - n1002ac4 + n1012ac4 + n1014ac4;
subject to dmbosord:
	0 <= pbosord0 + pbosord1 + pbosord2 + pbosord3 + pbosord4 + 302.0*-1 <= 61.0;
subject to dmboslga:
	0 <= pboslga0 + pboslga1 + pboslga2 + pboslga3 + 2352.0*-1 <= 471.0;
subject to dmboscle:
	0 <= pboscle0 + pboscle1 + pboscle2 + pboscle3 + 142.0*-1 <= 29.0;
subject to dmordbos:
	0 <= pordbos0 + pordbos1 + pordbos2 + pordbos3 + pordbos4 + 302.0*-1 <= 61.0;
subject to dmordlga:
	0 <= pordlga0 + pordlga1 + pordlga2 + pordlga3 + 515.0*-1 <= 103.0;
subject to dmordcle:
	0 <= pordcle0 + pordcle1 + pordcle2 + pordcle3 + 619.0*-1 <= 124.0;
subject to dmlgabos:
	0 <= plgabos0 + plgabos1 + plgabos2 + plgabos3 + plgabos4 + plgabos5 + 
	2743.0*-1 <= 549.0;
subject to dmlgaord:
	0 <= plgaord0 + plgaord1 + plgaord2 + plgaord3 + 712.0*-1 <= 143.0;
subject to dmlgacle:
	0 <= plgacle0 + plgacle1 + plgacle2 + 517.0*-1 <= 104.0;
subject to dmclebos:
	0 <= pclebos0 + pclebos1 + pclebos2 + pclebos3 + pclebos4 + 131.0*-1 <= 27.0;
subject to dmcleord:
	0 <= pcleord0 + pcleord1 + pcleord2 + pcleord3 + pcleord4 + 712.0*-1 <= 143.0;
subject to dmclelga:
	0 <= pclelga0 + pclelga1 + pclelga2 + pclelga3 + 409.0*-1 <= 82.0;
subject to msbosord:
	0 <= n1003ac1 + n1003ac2 + n1011ac1 + n1011ac2 + n1013ac2 + n1100ac2 + n1102ac2 
	+ 3.0*-1;
subject to msboslga:
	0 <= n1013ac2 + n1100ac2 + n1102ac2 + n1200ac2 + 7.0*-1;
subject to msboscle:
	0 <= n1005ac1 + n1005ac2 + n1011ac1 + n1011ac2 + n1100ac2 + n1102ac2 + -1;
subject to msordbos:
	0 <= n1006ac1 + n1006ac2 + n1008ac2 + n1010ac2 + n1100ac2 + n1102ac2 + 3.0*-1;
subject to msordlga:
	0 <= n1004ac2 + n1008ac2 + n1010ac2 + n1102ac2 + 4.0*-1;
subject to msordcle:
	0 <= n1002ac1 + n1002ac2 + n1006ac1 + n1006ac2 + n1010ac2 + n1100ac2 + 5.0*-1;
subject to mslgabos:
	0 <= n1008ac2 + n1010ac2 + n1014ac2 + n1100ac2 + n1102ac2 + n1201ac2 + 7.0*-1;
subject to mslgaord:
	0 <= n1013ac2 + n1017ac2 + n1019ac2 + n1100ac2 + 5.0*-1;
subject to mslgacle:
	0 <= n1015ac2 + n1019ac2 + n1100ac2 + 4.0*-1;
subject to msclebos:
	0 <= n1006ac1 + n1006ac2 + n1010ac2 + n1014ac2 + n1100ac2 + n1102ac2 + -1;
subject to mscleord:
	0 <= n1011ac1 + n1011ac2 + n1019ac2 + n1021ac1 + n1021ac2 + n1100ac2 + n1102ac2 
	+ 6.0*-1;
subject to msclelga:
	0 <= n1010ac2 + n1012ac2 + n1014ac2 + n1102ac2 + 3.0*-1;
subject to dcbosord:
	0 <= bbosord0 + bbosord1 + cbosord0 + cbosord1 + cbosord2 + cbosord3 + cbosord4 
	+ 12.0*-1 <= 12.0;
subject to dcboscle:
	0 <= bboscle0 + bboscle1 + cboscle0 + cboscle1 + cboscle2 + cboscle3 + 
	16.0*-1 <= 3.2;
subject to dcordbos:
	0 <= bordbos0 + cordbos0 + cordbos1 + cordbos2 + cordbos3 + cordbos4 + 
	24.0*-1 <= 4.8;
subject to dcordlga:
	0 <= cordlga0 + cordlga1 + cordlga2 + cordlga3 + 13.0*-1 <= 2.6;
subject to dclgaord:
	0 <= clgaord0 + clgaord1 + clgaord2 + clgaord3 + 45.0*-1 <= 9.0;
subject to dclgacle:
	0 <= clgacle0 + clgacle1 + clgacle2 + 16.0*-1 <= 3.2;
subject to dcclelga:
	0 <= cclelga0 + cclelga1 + cclelga2 + cclelga3 + 5.0*-1 <= 5.0;
subject to mcordbos:
	0 <= n1006ac3 + n1006ac4 + n1008ac4 + n1010ac4 + n1100ac4 + n1102ac4 + -1;
subject to mclgaord:
	0 <= n1013ac4 + n1017ac4 + n1019ac4 + n1100ac4 + 2.0*-1;

solve;
	display pbosord0;
	display pbosord1;
	display pbosord2;
	display pbosord3;
	display pbosord4;
	display pboslga0;
	display pboslga1;
	display pboslga2;
	display pboslga3;
	display pboscle0;
	display pboscle1;
	display pboscle2;
	display pboscle3;
	display pordbos0;
	display pordbos1;
	display pordbos2;
	display pordbos3;
	display pordbos4;
	display pordlga0;
	display pordlga1;
	display pordlga2;
	display pordlga3;
	display pordcle0;
	display pordcle1;
	display pordcle2;
	display pordcle3;
	display plgabos0;
	display plgabos1;
	display plgabos2;
	display plgabos3;
	display plgabos4;
	display plgabos5;
	display plgaord0;
	display plgaord1;
	display plgaord2;
	display plgaord3;
	display plgacle0;
	display plgacle1;
	display plgacle2;
	display pclebos0;
	display pclebos1;
	display pclebos2;
	display pclebos3;
	display pclebos4;
	display pcleord0;
	display pcleord1;
	display pcleord2;
	display pcleord3;
	display pcleord4;
	display pclelga0;
	display pclelga1;
	display pclelga2;
	display pclelga3;
	display bbosord0;
	display bbosord1;
	display bboscle0;
	display bboscle1;
	display bordbos0;
	display cbosord0;
	display cbosord1;
	display cbosord2;
	display cbosord3;
	display cbosord4;
	display cboscle0;
	display cboscle1;
	display cboscle2;
	display cboscle3;
	display cordbos0;
	display cordbos1;
	display cordbos2;
	display cordbos3;
	display cordbos4;
	display cordlga0;
	display cordlga1;
	display cordlga2;
	display cordlga3;
	display clgaord0;
	display clgaord1;
	display clgaord2;
	display clgaord3;
	display clgacle0;
	display clgacle1;
	display clgacle2;
	display cclelga0;
	display cclelga1;
	display cclelga2;
	display cclelga3;
	display grdtimo1;
	display grdtimn1;
	display grdtimo2;
	display grdtimn2;
	display grdtimn3;
	display grdtimn4;
	display n1003ac1;
	display n1003ac2;
	display n1003ac3;
	display n1003ac4;
	display n1005ac1;
	display n1005ac2;
	display n1005ac3;
	display n1005ac4;
	display n1011ac1;
	display n1011ac2;
	display n1011ac3;
	display n1011ac4;
	display n1013ac2;
	display n1013ac4;
	display n1015ac2;
	display n1015ac4;
	display n1017ac2;
	display n1017ac4;
	display n1019ac2;
	display n1019ac4;
	display n1021ac1;
	display n1021ac2;
	display n1021ac3;
	display n1021ac4;
	display n1002ac1;
	display n1002ac2;
	display n1002ac3;
	display n1002ac4;
	display n1004ac2;
	display n1004ac4;
	display n1006ac1;
	display n1006ac2;
	display n1006ac3;
	display n1006ac4;
	display n1008ac2;
	display n1008ac4;
	display n1010ac2;
	display n1010ac4;
	display n1012ac2;
	display n1012ac4;
	display n1014ac2;
	display n1014ac4;
	display n1100ac2;
	display n1100ac4;
	display n1102ac2;
	display n1102ac4;
	display n1200ac2;
	display n1200ac4;
	display n1201ac2;
	display n1201ac4;
display obj;
