(* ::Package:: *)

BeginPackage["Migdal`"];


Print["------------------------------------------
Migdal ionisation probabilities
P. Cox, M. Dolan, C. McCabe, H. Quiney (2022)
arXiv:2208.12222
------------------------------------------"];


$MigdalDirectory::usage="Location of Migdal ionisation tables.";

$alpha::usage="Fine-structure constant.";
$emin::usage="Lowest continuum electron energy in eV.";
$emax::usage="Highest continuum electron energy in eV.";
$vmin::usage="Lowest recoil velocity in units of c.";
$vmax::usage="Lowest recoil velocity in units of c.";

$elements::usage="List of implemented elements.";
$orbitals::usage="Dictionary of occupied orbitals for each atom.";

dpI1::usage="Differential probability for single ionisation without excitation.";
dpIE::usage="Differential probability for single ionisation with excitation.";
dpI21::usage="Differential probability for double ionisation with one electron above threshold.";
dpI2::usage="Differential probability for double ionisation with both electrons above threshold.";

p0::usage="Probability for no transition.";
pE::usage="Probability for excitation.";
pE1::usage="Probability for single excitation.";
pE2::usage="Probability for double excitation.";
pI1::usage="Integrated probability for single ionisation without excitation.";
pIE::usage="Integrated probability for single ionisation with excitation.";
pI21::usage="Integrated probability for double ionisation with one electron above threshold.";
pI2::usage="Integrated probability for double ionisation with both electrons above threshold.";
ptotal::usage="Total integrated probability.";

dpI1dipole::usage="Differential single ionisation probability using dipole approximation"
pI1dipole::usage"Integrated single ionisation probability using dipole approximation."


LoadTotalProbabilities::usage="LoadTotalProbabilities[element]
Load Migdal tables and define functions for total differential or integrated probabilities."

LoadProbabilities::usage="LoadProbabilities[element]
Load Migdal tables and define functions for differential or integrated probabilites."


Begin["`Private`"];


MigdalInterpolate::usage="MigdalInterpolate[filename]
Read Migdal table and interpolate.";


MigdalReset::usage="MigdalReset[]
Clear total probabilities";


MigdalResetOrbitals::usage="MigdalResetOrbitals[]
Clear orbital probabilities";


$MigdalDirectory=DirectoryName[$InputFileName];
$alpha = 1/137.0359895;
$emin = 1.0 10^-4; (* eV *)
$emax = 20.0; (* eV *)


$elements = {"He","C","F","Ne","Si","Ar","Ge","Kr","Xe"};

(* Occupied orbitals and binding energies (keV) *)
$orbitals = Association[{}];
AssociateTo[$orbitals,"He"->{{"1s",2.497980*^-02}}];
AssociateTo[$orbitals,"C"->{{"1s",3.083175*^-01},{"2s",1.921393*^-02},{"2p-",1.178626*^-02},{"2p",1.179453*^-02}}];
AssociateTo[$orbitals,"F"->{{"1s",7.186886*^-01},{"2s",4.288152*^-02},{"2p-",1.966527*^-02},{"2p",1.997926*^-02}}];
AssociateTo[$orbitals,"Ne"->{{"1s",8.930084*^-01},{"2s",5.267702*^-02},{"2p-",2.320668*^-02},{"2p",2.308252*^-02}}];
AssociateTo[$orbitals,"Si"->{{"1s",1.877873*^+00},{"2s",1.687017*^-01},{"2p-",1.167549*^-01},{"2p",1.159407*^-01},{"3s",1.500285*^-02},{"3p-",7.254801*^-03},{"3p",6.468341*^-03}}];
AssociateTo[$orbitals,"Ar"->{{"1s",3.241600*^+00},{"2s",3.377363*^-01},{"2p-",2.620990*^-01},{"2p",2.597887*^-01},{"3s",3.500978*^-02},{"3p-",1.620128*^-02},{"3p",1.599535*^-02}}];
AssociateTo[$orbitals,"Ge"->{{"1s",1.118596*^+01},{"2s",1.454912*^+00},{"2p-",1.288286*^+00},{"2p",1.256025*^+00},{"3s",2.019073*^-01},{"3p-",1.452096*^-01},{"3p",1.405931*^-01},{"3d-",4.426870*^-02},{"3d",4.359360*^-02},{"4s",1.566679*^-02},{"4p-",7.047754*^-03},{"4p",6.185060*^-03}}];
AssociateTo[$orbitals,"Kr"->{{"1s",1.441346*^+01},{"2s",1.961391*^+00},{"2p-",1.765333*^+00},{"2p",1.711030*^+00},{"3s",3.054330*^-01},{"3p-",2.345592*^-01},{"3p",2.262024*^-01},{"3d-",1.027950*^-01},{"3d",1.014111*^-01},{"4s",3.232023*^-02},{"4p-",1.473540*^-02},{"4p",1.399614*^-02}}];
AssociateTo[$orbitals,"Xe"->{{"1s",3.475594*^+01},{"2s",5.509354*^+00},{"2p-",5.161449*^+00},{"2p",4.835587*^+00},{"3s",1.170374*^+00},{"3p-",1.024780*^+00},{"3p",9.612494*^-01},{"3d-",7.081319*^-01},{"3d",6.948998*^-01},{"4s",2.293898*^-01},{"4p-",1.755814*^-01},{"4p",1.628001*^-01},{"4d-",7.377911*^-02},{"4d",7.166829*^-02},{"5s",2.748725*^-02},{"5p-",1.340357*^-02},{"5p",1.196770*^-02}}];


MigdalInterpolate::readerr="Error reading Migdal tables.";
Options[MigdalInterpolate] = {"EThreshold"->None, "Integrated"->False, "InterpolationOrder"->2, "Method"->"Spline", "PrecisionGoal"->6};
SyntaxInformation[MigdalInterpolate] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

MigdalInterpolate[filename_,type_,OptionsPattern[]]:= Module[{i,data,diff,diff2,EThreshold,int,int2,v0},
	If[OptionValue["EThreshold"] === None,
		EThreshold = $emin;
	, 
		EThreshold = OptionValue["EThreshold"];
		If[EThreshold < $emin, 
			EThreshold = $emin;
			OptionValue["EThreshold"] = None;
		];
	];

	data = Import[$MigdalDirectory<>filename,"Table"];
	If[data == $Failed,
		Message[MigdalInterpolate::readerr];
		Return[False];
	];
	
	Which[
		(* No transition *)
		type==0,
		data = {Log[#[[1]]],Log[#[[2]]]}& /@ data[[7;;]];
		int = Interpolation[data, Method->OptionValue["Method"], InterpolationOrder->1];
		Return[Exp[int[##]]&];
	,
		(* Excitation *)
		type==-1,
		data = {Log[#[[1]]],Log[#[[2]]]}& /@ data[[8;;]];
		int = Interpolation[data, Method->OptionValue["Method"], InterpolationOrder->OptionValue["InterpolationOrder"]];
		Return[Exp[int[##]]&];
	,
		(* Dipole single ionisation *)
		type==11,
		v0 = data[[11,2]];
		data = {Log[#[[1]]],Log[#[[3]]]}& /@ data[[11;;]];
		diff = Interpolation[data, Method->OptionValue["Method"], InterpolationOrder->OptionValue["InterpolationOrder"]];
		If[OptionValue["Integrated"],
			int = NIntegrate[Exp[lnEn] Exp[diff[lnEn]],{lnEn,Log[EThreshold],Log[$emax]},Method->"InterpolationPointsSubdivision",PrecisionGoal->OptionValue["PrecisionGoal"]];
			Return[int (Exp[#]/v0)^2&];
		,
			Return[(Exp[#2]/v0)^2 Exp[diff[#1]]&];
		];
	,
		(* Single ionisation or ionisation+excitation *)
		type==1 || type==-2,
		If[type==1, start=11, start=12];
		$vmin = data[[start,2]];
		$vmax = data[[-1,2]];
		If[OptionValue["Integrated"],
			data = {Log[#[[2]]],{Log[#[[1]]],#[[3]]}}& /@ data[[start;;]];
			data = Gather[data,#1[[1]]==#2[[1]]&];
			data = Table[{data[[i,1,1]],Interpolation[data[[i,All,2]], Method->OptionValue["Method"], InterpolationOrder->OptionValue["InterpolationOrder"]]},{i,1,Length[data]}];
			int = {#[[1]], Log[NIntegrate[Exp[lnEn] #[[2]][lnEn],{lnEn,Log[EThreshold],Log[$emax]},Method->"InterpolationPointsSubdivision",PrecisionGoal->OptionValue["PrecisionGoal"]]]}& /@ data;
			int = Interpolation[int, Method->OptionValue["Method"], InterpolationOrder->OptionValue["InterpolationOrder"]];
			Return[Exp[int[#]]&];
		,
			data = {{Log[#[[1]]],Log[#[[2]]]},Log[#[[3]]]}& /@ data[[start;;]];
			diff = Interpolation[data, Method->OptionValue["Method"], InterpolationOrder->OptionValue["InterpolationOrder"]];
			Return[Exp[diff[##]]&];
		];
	,
		(* Double ionisation *)
		type==2,
		If[OptionValue["Integrated"],
			(* Fully integrated probabilities *)
			data = {Log[#[[3]]],{{Log[#[[1]]],Log[#[[2]]]},#[[4]]}}& /@ data[[11;;]];
			data = Gather[data,#1[[1]]==#2[[1]]&];
			data = Table[{data[[i,1,1]],Interpolation[data[[i,All,2]], Method->OptionValue["Method"], InterpolationOrder->OptionValue["InterpolationOrder"]]},{i,1,Length[data]}];
			
			(* One electron above threshold *)
			If[!SameQ[OptionValue["EThreshold"], None],
				int = {#[[1]], Log[2*NIntegrate[Exp[lnEn1]Exp[lnEn2] #[[2]][lnEn1,lnEn2],{lnEn1,Log[$emin],Log[EThreshold]},{lnEn2,Log[EThreshold],Log[$emax]},Method->{"InterpolationPointsSubdivision",Method->"LocalAdaptive"},PrecisionGoal->OptionValue["PrecisionGoal"]]]}& /@ data;
				int = Interpolation[int, Method->OptionValue["Method"], InterpolationOrder->OptionValue["InterpolationOrder"]];
			];
		
			(* Two electrons above threshold *)			
			int2 = {#[[1]], Log[NIntegrate[Exp[lnEn1]Exp[lnEn2] #[[2]][lnEn1,lnEn2],{lnEn1,Log[EThreshold],Log[$emax]},{lnEn2,Log[EThreshold],Log[$emax]},Method->{"InterpolationPointsSubdivision",Method->"LocalAdaptive"},PrecisionGoal->OptionValue["PrecisionGoal"]]]}& /@ data;
			int2 = Interpolation[int2, Method->OptionValue["Method"], InterpolationOrder->OptionValue["InterpolationOrder"]];

			If[OptionValue["EThreshold"] === None,
				Return[{0&, Exp[int2[#]]&}];
				,
				Return[{Exp[int[#]]&, Exp[int2[#]]&}];
			];
		,
			(* Differential probabilities *)
			data = {{Log[#[[1]]],Log[#[[2]]],Log[#[[3]]]},Log[#[[4]]]}& /@ data[[11;;]];
			diff2 = Interpolation[data, Method->OptionValue["Method"], InterpolationOrder->OptionValue["InterpolationOrder"]];
			
			If[OptionValue["EThreshold"] === None,
				Return[{0&, Exp[diff2[##]]&}];
				,
				Return[{If[#1>Log[EThreshold],2*NIntegrate[Exp[lnEn2] Exp[diff2[#1,lnEn2,#2]],{lnEn2,Log[$emin],Log[EThreshold]},Method->{"InterpolationPointsSubdivision",Method->"LocalAdaptive"},PrecisionGoal->OptionValue["PrecisionGoal"]],0]&, Exp[diff2[##]]&}];
			];
		];
	];
];


Options[LoadProbabilities]=Flatten[{Options[MigdalInterpolate],{"DarkMatter"->False, "Dipole"->False, "Double"->False, "Inclusive"->False, "VelocityGrid"->"linear"}}];
SyntaxInformation[LoadProbabilities] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

LoadProbabilities[element_,opts:OptionsPattern[]]:=Module[{i, j, pair, ppI1, ppI2},
	MigdalReset[Dipole->OptionValue["Dipole"], Integrated->OptionValue["Integrated"]];

	If[!MemberQ[$elements,element],
		Print["Tables for "<>element<>" are not available."];
		Return[1];
	];

	If[OptionValue["Inclusive"] && OptionValue["Dipole"],
		Print["Dipole and inclusive are incompatible options."];
		Return[];
	]; 
			
	If[OptionValue["Inclusive"] && (OptionValue["EThreshold"]===None || OptionValue["EThreshold"] < 0.5),
		Print["EThreshold must be at least 0.5 keV for semi-inclusive probabilities."];
		Return[];
	]; 

	If[OptionValue["Integrated"],
		PrintTemporary[" Loading integrated probabilities for "<>element<>"..."];
	,
		PrintTemporary[" Loading differential probabilities for "<>element<>"..."];
	];
	
	(* Dipole single ionisation *)
	If[OptionValue["Dipole"],
		ppI1 = Association[#[[1]]->MigdalInterpolate[element<>"/single-ionisation/"<>element<>"_"<>#[[1]]<>"_dipole.txt",11,FilterRules[{opts}, Options[MigdalInterpolate]]]& /@ $orbitals[[element]]];
		If[OptionValue["Integrated"], 
			pI1dipole = ppI1;
			Print["Loading integrated probabilities for "<>element<>"...done."];
		, 
			dpI1dipole = ppI1;
			Print["Loading differential probabilities for "<>element<>"...done."];
		];
		Return[];
	];
	
	(* Semi-inclusive ionisation *)
	If[OptionValue["Inclusive"], 
		ppI1 = Association[#[[1]]->MigdalInterpolate[element<>"/semi-inclusive/"<>element<>"_"<>#[[1]]<>"_semi-inclusive.txt",1,FilterRules[{opts}, Options[MigdalInterpolate]]]& /@ $orbitals[[element]]];
		If[OptionValue["Integrated"], 
			pI1 = ppI1;
			Print["Loading integrated probabilities for "<>element<>"...done."];
		, 
			dpI1 = ppI1;
			Print["Loading differential probabilities for "<>element<>"...done."];
		];
		Return[];
	];
	
	(* Single ionisation *)
	If[OptionValue["DarkMatter"] || ToLowerCase[OptionValue["VelocityGrid"]]=="log",
		(* Logarithmic velocity grid *)
		ppI1 = Association[#[[1]]->MigdalInterpolate[element<>"/single-ionisation/"<>element<>"_"<>#[[1]]<>"_DM.txt",1,FilterRules[{opts}, Options[MigdalInterpolate]]]& /@ $orbitals[[element]]];
	,
		(* Linear velocity grid *)
		ppI1 = Association[#[[1]]->MigdalInterpolate[element<>"/single-ionisation/"<>element<>"_"<>#[[1]]<>".txt",1,FilterRules[{opts}, Options[MigdalInterpolate]]]& /@ $orbitals[[element]]];
	];
	If[OptionValue["Integrated"], pI1 = ppI1, dpI1 = ppI1];
	
	(* Double ionisation *)
	If[OptionValue["Double"] && !OptionValue["Inclusive"],
		For[i=1, i<=Length[$orbitals[[element]]], i++,
			For[j=i, j<=Length[$orbitals[[element]]], j++,
				pair = $orbitals[[element]][[i]][[1]] <> $orbitals[[element]][[j]][[1]];
				If[element=="C" && pair == "2p-2p", Continue[]];
				
				ppI2 = MigdalInterpolate[element<>"/double-ionisation/"<>element<>"_"<>pair<>".txt",2,FilterRules[{opts}, Options[MigdalInterpolate]]];
				
				If[OptionValue["Integrated"],
					AssociateTo[pI21, pair->ppI2[[1]]]; 
					AssociateTo[pI2, pair->ppI2[[2]]];
				,
					AssociateTo[dpI21, pair->ppI2[[1]]];
					AssociateTo[dpI2, pair->ppI2[[2]]];
				];
			];
		];	
	];
	
	If[OptionValue["Integrated"],
		Print["Loading integrated probabilities for "<>element<>"...done."];
	,
		Print["Loading differential probabilities for "<>element<>"...done."];
	];
];


Options[LoadTotalProbabilities]=Flatten[{Options[MigdalInterpolate],{"DarkMatter"->False, "Dipole"->False, "Double"->False, "Excitations"->False, "Inclusive"->False, "VelocityGrid"->"linear"}}];
SyntaxInformation[LoadTotalProbabilities] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

LoadTotalProbabilities[element_,opts:OptionsPattern[]]:= Module[{ppI1,ppIE,ppI2},
	MigdalResetTotal[Dipole->OptionValue["Dipole"], Integrated->OptionValue["Integrated"]];
	
	If[!MemberQ[$elements,element],
		Print["Tables for "<>element<>" are not available."];
		Return[1];
	];
	
	If[OptionValue["Inclusive"] && OptionValue["Dipole"],
		Print["Dipole and inclusive are incompatible options."];
		Return[];
	]; 
	
	If[OptionValue["Inclusive"] && (OptionValue["EThreshold"]===None || OptionValue["EThreshold"] < 0.5),
		Print["EThreshold must be at least 0.5 keV for semi-inclusive probabilities."];
		Return[];
	]; 
	
	If[OptionValue["Integrated"],
		PrintTemporary[" Evaluating integrated probabilities for "<>element<>"..."];
	,
		PrintTemporary[" Loading differential probabilities for "<>element<>"..."];
	];
	
	(* Dipole single ionisation *)
	If[OptionValue["Dipole"],
		ppI1 = MigdalInterpolate[element<>"/single-ionisation/"<>element<>"_single-ionisation_dipole"<>".txt",11,FilterRules[{opts}, Options[MigdalInterpolate]]];
		If[OptionValue["Integrated"], 
			pI1dipole = ppI1;
			Print["Evaluating integrated probabilities for "<>element<>"...done."];
		, 
			dpI1dipole = ppI1;
			Print["Evaluating differential probabilities for "<>element<>"...done."];
		];
		Return[];
	];
	
	(* No transition *)
	p0 = MigdalInterpolate[element<>"/"<>element<>"_no-transition.txt",0];

	(* Semi-inclusive ionisation *)
	If[OptionValue["Inclusive"], 
		ppI1 = MigdalInterpolate[element<>"/semi-inclusive/"<>element<>"_semi-inclusive.txt",1,FilterRules[{opts}, Options[MigdalInterpolate]]];
		If[OptionValue["Integrated"], 
			pI1 = ppI1;
			Print["Evaluating integrated probabilities for "<>element<>"...done."];
		,
			dpI1 = ppI1;
			Print["Evaluating differential probabilities for "<>element<>"...done."];
		];
		Return[];
	];
		
	(* Single ionisation *)
	If[OptionValue["DarkMatter"] || OptionValue["VelocityGrid"]=="log",
		(* Logarithmic velocity grid *)
		ppI1 = MigdalInterpolate[element<>"/single-ionisation/"<>element<>"_single-ionisation_DM.txt",1,FilterRules[{opts}, Options[MigdalInterpolate]]];
	,
		(* Linear velocity grid *)
		ppI1 = MigdalInterpolate[element<>"/single-ionisation/"<>element<>"_single-ionisation.txt",1,FilterRules[{opts}, Options[MigdalInterpolate]]];
	];
	If[OptionValue["Integrated"], pI1 = ppI1, dpI1 = ppI1];
	
	(* Excitations *)
	If[OptionValue["Excitations"],
		
		(* Single excitation *)
		pE1 = MigdalInterpolate[element<>"/excitations/"<>element<>"_single-excitation.txt",-1];
		pE = pE1;
			
		If[OptionValue["Double"],
			(* Double excitation *)
			pE2 = MigdalInterpolate[element<>"/excitations/"<>element<>"_double-excitation.txt",-1];
			pE = pE1[#]+pE2[#]&;
				
			(* Single ionisation + excitation *)
			ppIE = MigdalInterpolate[element<>"/excitations/"<>element<>"_excitation+ionisation.txt",-2,FilterRules[{opts}, Options[MigdalInterpolate]]];
			If[OptionValue["Integrated"], pIE = ppIE, dpIE = ppIE];
		];
	];
	
	(* Double ionisation *)
	If[OptionValue["Double"],
		ppI2 = MigdalInterpolate[element<>"/double-ionisation/"<>element<>"_double-ionisation.txt",2,FilterRules[{opts}, Options[MigdalInterpolate]]];
		If[OptionValue["Integrated"], 
			pI21 = ppI2[[1]];
			pI2 = ppI2[[2]];
		, 
			dpI21 = ppI2[[1]];
			dpI2 = ppI2[[2]];
		];
	];
	
	(* Total probabilities *)
	If[OptionValue["Integrated"], ptotal[lnv_]:= p0[lnv] + pE[lnv] + pI1[lnv] + pIE[lnv] + pI2[lnv];];
	
	If[OptionValue["Integrated"],
		Print["Evaluating integrated probabilities for "<>element<>"...done."];
	,
		Print["Evaluating differential probabilities for "<>element<>"...done."];
	];
];


MigdalReset::usage="MigdalReset[]
Set probabilities to zero."
Options[MigdalReset]= {"Dipole"->False, "Integrated"->False};
SyntaxInformation[MigdalReset] = {"ArgumentsPattern" -> {OptionsPattern[]}};

MigdalReset[OptionsPattern[]]:=Module[{},

	If[OptionValue["Dipole"],
		If[OptionValue["Integrated"],
			Clear[pI1dipole];
			pI1dipole = Association[{}];
		,
			Clear[dpI1dipole];
			dpI1dipole = Association[{}];
		];
		Return[];
	];

	If[OptionValue["Integrated"],
		Clear[pI1,pI21,pI2,ptotal];
		pI1 = Association[{}];
		pI21 = Association[{}];
		pI2 = Association[{}];
	,
		Clear[dpI1,dpI21,dpI2];
		dpI1 = Association[{}];
		dpI21 = Association[{}];
		dpI2 = Association[{}];
	];
];


MigdalResetTotal::usage="MigdalResetTotal[]
Set probabilities to zero."
Options[MigdalResetTotal]= {"Dipole"->False,"Integrated"->False};
SyntaxInformation[MigdalResetTotal] = {"ArgumentsPattern" -> {OptionsPattern[]}};

MigdalResetTotal[OptionsPattern[]]:=Module[{},

	If[OptionValue["Dipole"],
		If[OptionValue["Integrated"],
			Clear[pI1dipole];
			pI1dipole[lnv_]:= 0;
		,
			Clear[dpI1dipole];
			dpI1dipole[lnv_,lnEn_]:= 0;
		];
		Return[];
	];

	Clear[p0,pE,pE1,pE2];
	p0[lnv_]:= 0;
	pE1[lnv_]:= 0;
	pE2[lnv_]:= 0;
	pE[lnv_]:= 0;

	If[OptionValue["Integrated"],
		Clear[pI1,pIE,pI21,pI2,ptotal];
		pI1[lnv_]:= 0;
		pIE[lnv_]:= 0;
		pI21[lnv_]:= 0;
		pI2[lnv_]:= 0;
		ptotal[lnv_]:= 0;
	,
		Clear[dpI1,dpIE,dpI21,dpI2];
		dpI1[lnv_,lnEn_]:= 0;
		dpIE[lnv_,lnEn_]:= 0;
		dpI21[lnv_,lnEn_]:= 0;
		dpI2[lnv_,lnEn1_,lnEn2_]:= 0;
	];
];


End[];


EndPackage[];
