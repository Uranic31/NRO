(* ::Package:: *)

priblizekPi[stTock_]:=Module[{tocke,krogNot,kroznica},
tocke=tockeN[stTock];
krogNot=Select[tocke,Norm[#] <= 1 &];
kroznica=Select[tocke,krog];
4 (Length[krogNot] + Length[kroznica])/Length[tocke]]
