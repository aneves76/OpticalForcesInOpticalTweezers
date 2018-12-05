(* ::Package:: *)

(* ::Program:: *)
BeginPackage["MieCoefficients`"]


ABMieCoefficients::usage="ABMieCoefficients[n,x,m] returns the standard Mie Coefficients, in the form {A,B} where n is the order, x the size factor (defined as x=k a, where k is the wavenumber in the surrounding medium, and a is the sphere radius), m is the relative refractive index (m=m_particle/m_surroundingmedium)."

Begin["`Private`"]

ABMieCoefficients[n_,x_,m_]=Module[
{\[Psi],\[Xi],D\[Psi],D\[Xi]},
\[Psi][n_,x_]=x SphericalBesselJ[n,x];
\[Xi][n_,x_]=x SphericalHankelH1[n,x];
D\[Psi][n_,x_]=1/2 (x SphericalBesselJ[-1+n,x]+SphericalBesselJ[n,x]-x SphericalBesselJ[1+n,x]);
D\[Xi][n_,x_]=1/2 (x SphericalHankelH1[-1+n,x]+SphericalHankelH1[n,x]-x SphericalHankelH1[1+n,x]);
{(m \[Psi][n,m x] D\[Psi][n,x]-\[Psi][n,x] D\[Psi][n,m x])/(\[Xi][n,x] D\[Psi][n,m x]-m \[Psi][n,m x] D\[Xi][n,x]),(m D\[Psi][n,m x] \[Psi][n,x]-D\[Psi][n,x] \[Psi][n,m x])/(D\[Xi][n,x] \[Psi][n,m x]-m D\[Psi][n,m x] \[Xi][n,x])}
];

End[]

EndPackage[]
