(* ::Package:: *)

(* ::Program:: *)
(**)


BeginPackage["IncidentBeams`"]


GTMpw::usage="GTMpw[n,m] returns the transverse magnetic beam shape coeficient for a plane wave."
GTEpw::usage="GTEpw[n,m] returns the transverse eletric beam shape coeficient for a plane wave."
BSCfocGaussian::usage="BSCfocGaussian[n,m,k\[Rho]0,\[Phi]0,kz0,px,py,pre,\[Alpha]NA,f\[Omega]] returns the beam sape coeficient for a highly focused Gaussina beam, located at (k\[Rho]0,\[Phi]0,kz0) from the origin, where k is the wavenumber in the medium. The Jones vector correspontds to (px,py), the focusing system has a maximum numerical aperture of \[Alpha]NA, focal length f and beam waist at the back objective aperture of \[Omega]. f\[Omega] correspond to f/\[Omega], and pre to kf Exp[\[ImaginaryI]kf]\!\(\*SqrtBox[FractionBox[\(n_b\), \(n_a\)]]\), where n_b it the refractive index before and n_a after the focusing lens."

Begin["`Private`"]

(*Define the BSC for x-polarized plane wave*)
GTMpw[n_Integer?Positive,m_Integer]=-(KroneckerDelta[m,1]-KroneckerDelta[m,-1])Sqrt[\[Pi](2 n+1)] I^(n+1);
GTEpw[n_Integer?Positive,m_Integer]=-I(KroneckerDelta[m,1]+KroneckerDelta[m,-1])Sqrt[\[Pi](2 n+1)] I^(n+1);

(*Define the BSC for highly focused Gaussian beam*)
(*Returns {GTM,GTE} which corresponds to Transverse Magnetic and Eletric BSC for Jones vector {p_x,p_y}*)
(*pre corresponds to kf Exp[\[ImaginaryI]kf]Sqrt[n_b/n_a]*)
(*f\[Omega] corresponds to f/\[Omega]*)
BSCfocGaussian[n_Integer?Positive,m_Integer,k\[Rho]0_?Positive,\[Phi]0_,kz0_,px_,py_,pre_,\[Alpha]NA_,f\[Omega]_] /; Abs[m] <= n:=Module[
{G1,G2},
G1=-Sqrt[(4\[Pi](2n+1))/(n(n+1))]I^(n-m+1) pre NIntegrate[Sqrt[(n-m)!/(n+m)!]Sin[\[Alpha]]Sqrt[Cos[\[Alpha]]]Exp[-f\[Omega]^2 Sin[\[Alpha]]^2]Exp[-I kz0 Cos[\[Alpha]]]((m BesselJ[m,k\[Rho]0 Sin[\[Alpha]]])/(k\[Rho]0 Sin[\[Alpha]]) (m LegendreP[n,m,Cos[\[Alpha]]])/Sin[\[Alpha]]-1/2 Csc[\[Alpha]] ((1+n) Cos[\[Alpha]] LegendreP[n,m,Cos[\[Alpha]]]+(-1+m-n) LegendreP[1+n,m,Cos[\[Alpha]]]) (BesselJ[-1+m,k\[Rho]0 Sin[\[Alpha]]]-BesselJ[1+m,k\[Rho]0 Sin[\[Alpha]]])),{\[Alpha],0,\[Alpha]NA},WorkingPrecision->25,MaxRecursion->100,PrecisionGoal->9,AccuracyGoal->8,Method->{"OscillatorySelection","SymbolicProcessing"->0}];
G2=-Sqrt[(4\[Pi](2n+1))/(n(n+1))]I^(n-m+1) pre NIntegrate[Sqrt[(n-m)!/(n+m)!]Sin[\[Alpha]]Sqrt[Cos[\[Alpha]]]Exp[-f\[Omega]^2 Sin[\[Alpha]]^2]Exp[-I kz0 Cos[\[Alpha]]](-Csc[\[Alpha]] ((1+n) Cos[\[Alpha]] LegendreP[n,m,Cos[\[Alpha]]]+(-1+m-n) LegendreP[1+n,m,Cos[\[Alpha]]]) (m BesselJ[m,k\[Rho]0 Sin[\[Alpha]]])/(k\[Rho]0 Sin[\[Alpha]])+1/2 (BesselJ[-1+m,k\[Rho]0 Sin[\[Alpha]]]-BesselJ[1+m,k\[Rho]0 Sin[\[Alpha]]]) (m LegendreP[n,m,Cos[\[Alpha]]])/Sin[\[Alpha]]),{\[Alpha],0,\[Alpha]NA},WorkingPrecision->25,MaxRecursion->100,PrecisionGoal->9,AccuracyGoal->8,Method->{"OscillatorySelection","SymbolicProcessing"->0}];
{G1(px Cos[\[Phi]0]+py Sin[\[Phi]0])+ G2(-px Sin[\[Phi]0]+py Cos[\[Phi]0]),G1(py Cos[\[Phi]0]-px Sin[\[Phi]0])+ G2(-py Sin[\[Phi]0]-px Cos[\[Phi]0])}
];

BSCfocGaussian[n_Integer?Positive,m_Integer,0,\[Phi]0_,kz0_,px_,py_,pre_,\[Alpha]NA_,f\[Omega]_] /; Abs[m] <= n:=Module[
{G},
G=Piecewise[{{(pre I^(n) Sqrt[\[Pi](2 n+1)])/(n(n+1)) NIntegrate[Sqrt[Cos[\[Alpha]]]Exp[-f\[Omega]^2 Sin[\[Alpha]]^2]Exp[-I kz0 Cos[\[Alpha]]](-(-1+(1+n) Cos[\[Alpha]]) LegendreP[n,1,Cos[\[Alpha]]]+n LegendreP[1+n,1,Cos[\[Alpha]]]),{\[Alpha],0,\[Alpha]NA},WorkingPrecision->25,MaxRecursion->100,PrecisionGoal->9,AccuracyGoal->8,Method->{"OscillatorySelection","SymbolicProcessing"->0}],m==1},{(pre I^(n) Sqrt[\[Pi](2 n+1)])/(n(n+1)) NIntegrate[Sqrt[Cos[\[Alpha]]]Exp[-f\[Omega]^2 Sin[\[Alpha]]^2]Exp[-I kz0 Cos[\[Alpha]]](-(-1+(1+n) Cos[\[Alpha]]) LegendreP[n,1,Cos[\[Alpha]]]+n LegendreP[1+n,1,Cos[\[Alpha]]]),{\[Alpha],0,\[Alpha]NA},WorkingPrecision->25,MaxRecursion->100,PrecisionGoal->9,AccuracyGoal->8,Method->{"OscillatorySelection","SymbolicProcessing"->0}],m==-1}},0];
G{((KroneckerDelta[m,-1] -KroneckerDelta[m,1]) px + I py),(-I(-KroneckerDelta[m,-1]+KroneckerDelta[m,1]) px - py)}
];

End[]

EndPackage[]
