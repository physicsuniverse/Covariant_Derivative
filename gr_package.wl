(* ::Package:: *)

(* :Title: General Relativity Package *)
(* :Context: GRPackage` *)
(* :Author: Peng Liu *)
(* :Date: June 12, 2025 *)
(* :Version: 1.0 *)
(* :Mathematica Version: 12.0+ *)
(* :Description: 
   A Mathematica package for symbolic general relativity calculations.
   Provides functions for computing Christoffel symbols, Riemann tensor,
   Ricci tensor, Ricci scalar, Einstein tensor, Weyl tensor, and 
   covariant derivatives using elegant symbolic operations.
*)

BeginPackage["GRPackage`"];

(* Usage messages *)
ChristoffelSymbol::usage = "ChristoffelSymbol[g, xx] computes the Christoffel symbols from metric tensor g and coordinates xx.";
RiemannTensor::usage = "RiemannTensor[g, xx] computes the Riemann curvature tensor using Wald's convention.";
RicciTensor::usage = "RicciTensor[g, xx] computes the Ricci tensor.";
RicciScalar::usage = "RicciScalar[g, xx] computes the Ricci scalar.";
EinsteinTensor::usage = "EinsteinTensor[g, xx] computes the Einstein tensor.";
WeylTensor::usage = "WeylTensor[g, xx] computes the Weyl conformal tensor.";
CovariantDerivative::usage = "CovariantDerivative[a, ten, g, xx] computes the covariant derivative of tensor ten with respect to coordinate a, using metric g and coordinates xx. Upper indices are denoted with minus sign (e.g., -a), following xAct notation.";
permuteTranpose::usage = "permuteTranpose[t, range, i, or] is a helper function for tensor permutations.";

Begin["`Private`"];

Clear[ChristoffelSymbol,RiemannTensor,RicciTensor,RicciScalar,EinsteinTensor,WeylTensor];

ChristoffelSymbol[g_,xx_]:=Block[
	{
		ig,gd
		},
	ig = Inverse[g];
	gd = D[g,{xx}];
	ig . (gd+Transpose[gd,{3,1,2}]-Transpose[gd,{2,3,1}])/2
];

RiemannTensor[g_,xx_]:=Block[
	{
		ig, chr, chrd, chrdif
	},

	chr    = ChristoffelSymbol[g, xx];
	chrd   = D[chr, {xx}];
	chrdif = chrd - chr . chr;
	Transpose[chrdif, {4, 1, 3, 2}]-Transpose[chrdif, {4, 2, 3, 1}]
];

(* RiemannTensor here adopts convention of Wald's <General Relativity> *)
RicciTensor[g_,xx_]:=Block[
	{
		rie
	},
	rie = RiemannTensor[g,xx];
	Tr[Transpose[rie,{3,2,4,1}],Plus,2]
];

RicciScalar[g_,xx_]:=Block[
	{ric},
	ric = RicciTensor[g,xx];
	ric . Inverse[g]//Tr
];

EinsteinTensor[g_,xx_]:=Block[
	{
		rict,rics
	},
	RicciTensor[g,xx]-RicciScalar[g,xx]g/2
];

WeylTensor[g_,xx_]:=Block[
	{
		rie,ric,rics,n,granti,gganti,part2,part3
	},

	n               = Length[xx];
	{rie,ric,rics}  = Through[{RiemannTensor,RicciTensor,RicciScalar}[g,xx]];
	{gganti,granti} = (#-Transpose[#,{1,3,2,4}])/2& /@ (TensorProduct[g,#]& /@ {g,ric});
	part2           = Transpose[granti,{1,3,4,2}]-Transpose[granti,{2,3,4,1}];
	part3           = rics*Transpose[gganti,{1,3,4,2}];
	rie . g-(2part2)/(n-2)+(2part3)/((n-1)(n-2))
];

permuteTranpose[t_,range_,i_,or_]:=Transpose[
	t,
	Permute[Range[range],Cycles[{Range[i]//If[or===0,Reverse,Identity]}]]
];

ClearAll[CovariantDerivative];
SetAttributes[CovariantDerivative,HoldAll];
CovariantDerivative[a_,ten_,g_,xx_]:=Block[
	{
		chr, tenhold, tenlist, nt, posi, part, t1, res, alist, tensorindices
	},
	
	chr     = ChristoffelSymbol[g, xx];
	tenlist = List @@ ten;
	tensorindices = StringSplit[
		StringTake[ToString[ReleaseHold[Hold@@@Hold[ten]]], {6, -2}], ", "
		] /. str_String /; StringTake[str, 1]==="-" :> -StringTake[str, {2,-1}];
	
	t1      = ten[[0]];
	nt      = Length@tenlist;

	posi[1] = Position[tenlist, _Times]//Flatten;
	posi[2] = Complement[Range[nt], posi[1]];

	part[0] = Transpose[D[t1,{xx}], RotateLeft[Range[nt+1]]];
	part[1] = Sum[
		permuteTranpose[chr . permuteTranpose[t1, nt, posi[1][[i]], 0], nt+1, posi[1][[i]]+1, 1],
		{i, Length@posi[1]}
	];
	part[2] = -Sum[
		permuteTranpose[Transpose[chr, {3, 1, 2}] . permuteTranpose[t1, nt, posi[2][[i]], 0], nt+1, posi[2][[i]]+1, 1],
		{i,Length@posi[2]}
	];

	res     = If[
		Head[a] === Times,
		Inverse[g] . Total[part /@ Range[0,2]],
		Total[part /@ Range[0,2]]
	];

	alist   = Flatten[{a,tenlist}] /. -aa_ :> aa;

	If[
		MemberQ[alist[[2;;-1]],alist[[1]]],
		TensorContract[res,Flatten[Position[alist,alist[[1]]],{2}]],
		res
	]
];

End[];
EndPackage[];