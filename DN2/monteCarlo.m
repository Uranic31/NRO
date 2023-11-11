(* ::Package:: *)

(* ::Input:: *)
(*tockeN[stTock_]:=Table[{RandomReal[],RandomReal[]},{stTock}]*)


(* ::Input:: *)
(*krog[{x_,y_}]:=x^2+y^2==1*)


(* ::Input:: *)
(*Get["C:\\Users\\lukau\\Desktop\\DN2\\priblizekPi.m"]*)


(* ::Input:: *)
(*stTock=1000;*)


(* ::Input:: *)
(*priblizek=priblizekPi[stTock];*)
(*priblizek=N[priblizek]*)


(* ::Input:: *)
(*tocke=tockeN[stTock];*)
(*tockeKrog=Select[tocke,Norm[#]<=1&];*)
(*kroznica=Select[tocke,krog];*)
(*zunaj=Complement[tocke,Join[tockeKrog,kroznica]];*)


(* ::Input:: *)
(*p1=ListPlot[tockeKrog,PlotStyle->Blue,AspectRatio->1,PlotRange->{{0,1},{0,1}},PlotLegends->{"To\[CHacek]ke v krogu"}];*)
(*p2=ListPlot[zunaj,PlotStyle->Red,PlotLegends->{"To\[CHacek]ke zunaj kroga"}];*)
(*p3=Graphics[{EdgeForm[Black],FaceForm[None],Circle[{0,0},1]},Axes->True,AspectRatio->1,PlotRange->{{0,1},{0,1}}];*)


(* ::Input:: *)
(*Show[p1,p2,p3,PlotRange->{{0,1},{0,1}}]*)
