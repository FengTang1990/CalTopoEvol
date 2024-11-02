(* ::Package:: *)

BeginPackage["CalTopoEvol`"];


dir::usage="specify the directory of source data: CalTopoEvol.mx."
indhspirreps::usage="This relates irreps of indhsp sym. data from differeint settings (however, here, it is not used)."
SBPattern::usage="symmetry-breaking patterns (0,1,2,etc.) for the 1651 MSGs."
allhspsinlowmsg::usage=" "
allhslhsplhspsinlowmsg::usage=" "
hsprelatedV2::usage=" "
allAIB::usage=" "
allhslbandnode::usage="not used"
allhspbandnode::usage="not used"
CalTopoSubGrp::usage=" "
CalTopoEvol::usage=" "
SIValues::usage=" "
ShowBC::usage=" "
allhslhsplhsps::usage="not used"
allimsg2::usage=" "
indhspname::usage=" "
resultfrom1toi::usage=" "
submpg::usage=" "
mpgmsg::usage=" "
allmpg::usage=" "
lowmsg::usage=" "
allmsg3::usage=" "
hsprelatedV2inlowmsgCor::usage=" "


Begin["`Private`"];


dir=Input["Please input the directory of CalTopoEvol.mx:"];
{indhspirreps,SBPattern,allhspsinlowmsg,allhspbandnode,hsprelatedV2,allAIB,
allhslbandnode,allhslhsplhspsinlowmsg,allhslhsplhsps,allimsg2,indhspname,resultfrom1toi,submpg,mpgmsg,allmpg,lowmsg
,allmsg3,hsprelatedV2inlowmsgCor}=Import[dir<>"\\CalTopoEvol.mx"];


CalTopoSubGrp[imsg_,symbreakpattern_,n_,so_]:=Module[{s1,s2,AIB,q,ilow,case,nu,i,ii,i1,i2,i3,i4,i5,deg,hsp,hsl,hspl,s3,nn,n3,cor,\[Delta]n,nnn,d3,jjj,tofill,permu,tmpnu,nnn0,nnn00,NN,NNA,tuples,bnhsp,bnhsl,bnhspl,tmptuples,tuplesA,tuplesB,RESULT,hsl2,hspl2,n1,n2,nn1,nn2,hsp0,hsp00,\[Delta]nn12,n7,n8,i33,dBS,XBS,basis,SI,iset,hsp2,
tmphsp1,tohsp1,f1025,tmp1101,tmp1101A,tmp1101B,tmp1101C},
iset=1;
s1=s2=0;
ilow=Position[SBPattern[[imsg]],symbreakpattern][[1,1]];
hsp=Table[allhspsinlowmsg[[imsg,iset,i1,i2,so+1,ilow]],{i1,Length[allhspsinlowmsg[[imsg,iset]]]},{i2,Length[allhspsinlowmsg[[imsg,iset,i1]]]}];
If[imsg>7,
hsl=Table[allhslhsplhspsinlowmsg[[imsg-7,iset,1,i1,i2,i3,ilow,so+1]],{i1,Length[allhslhsplhspsinlowmsg[[imsg-7,iset,1]]]},
{i2,Length[allhslhsplhspsinlowmsg[[imsg-7,iset,1,i1]]]},
{i3,Length[allhslhsplhspsinlowmsg[[imsg-7,iset,1,i1,i2]]]}];
hspl=Table[allhslhsplhspsinlowmsg[[imsg-7,iset,2,i1,i2,i3,ilow,so+1]],{i1,Length[allhslhsplhspsinlowmsg[[imsg-7,iset,2]]]},
{i2,Length[allhslhsplhspsinlowmsg[[imsg-7,iset,2,i1]]]},
{i3,Length[allhslhsplhspsinlowmsg[[imsg-7,iset,2,i1,i2]]]}],
hsl=hspl={}];
hsp0=Table[hsp[[i1,i2,1]],{i1,Length[hsp]},{i2,Length[hsp[[i1]]]}];
hsp00=Table[Last[hsp[[i1,i2]]],{i1,Length[hsp]},{i2,Length[hsp[[i1]]]}];
hsp2=hsprelatedV2[[imsg,iset]];
If[imsg>7,
{hsl2,hspl2}=allhslhsplhsps[[imsg-7,iset,1]];
,hsl2=hspl2={}];

deg=Table[hsp[[i1,1,{1,2}]],{i1,Length[hsp]}];
s1++;If[Length[deg]==Length[n],s2++];
nu={};
For[i=1,i<=Length[deg],i++,
s1++;
If[deg[[i,1]]==n[[i,1]],s2++];
s1++;If[n[[i,2]]\[Transpose][[2]]==deg[[i,2]],s2++];
AppendTo[nu,Sum[deg[[i,2,i1]]n[[i,2,i1,3]],{i1,Length[deg[[i,2]]]}]]];
nu=FullSimplify[nu];
s1++;If[Length[DeleteDuplicates[nu]]==1,s2++];
(*all hsp/hsl/hspl find band nodes*)
NN={};
For[i1=1,i1<=Length[hsp],i1++,NNA={};
For[i2=1,i2<=Length[hsp[[i1]]],i2++,
nn=hsp[[i1,i2,3]]\[Transpose] . n[[i1,2]]\[Transpose][[3]];
nn=FullSimplify[nn];s3={};
For[i3=1,i3<=Length[nn],i3++,
If[IntegerQ[nn[[i3]]],Continue[]];
AppendTo[s3,i3]];
s1++;If[Length[s3]<=1,s2++];
n3=s3;
cor=hsp[[i1,i2,5]];
If[Length[n3]==1,
\[Delta]n=nn[[n3[[1]]]]-Floor[nn[[n3[[1]]]]];
d3=hsp[[i1,i2,2,n3[[1]]]];
\[Delta]n=FullSimplify[\[Delta]n d3];
nn[[n3[[1]]]]=Floor[nn[[n3[[1]]]]];
nnn0=FullSimplify[cor\[Transpose] . nn];
jjj=cor[[n3[[1]]]];
tofill={};
For[i4=1,i4<=Length[jjj],i4++,
For[i5=1,i5<=jjj[[i4]],i5++,
AppendTo[tofill,{i4,hsp[[i1,i2,4,i4]]}]]];
permu=Permutations[Range[Length[tofill]]];
nnn={};
For[i4=1,i4<=Length[permu],i4++,
tmpnu=0;nnn00=nnn0;
For[i5=1,i5<=Length[permu[[i4]]],i5++,
tmpnu+=tofill[[permu[[i4,i5]],2]];
If[tmpnu>=\[Delta]n,Break[]];
nnn00[[tofill[[permu[[i4,i5]],1]]]]=nnn00[[tofill[[permu[[i4,i5]],1]]]]+1];
s1++;If[i5<=Length[permu[[i4]]],s2++];
nnn00[[tofill[[permu[[i4,i5]],1]]]]=nnn00[[tofill[[permu[[i4,i5]],1]]]]+(tofill[[permu[[i4,i5]],2]]-(tmpnu-\[Delta]n))/tofill[[permu[[i4,i5]],2]];
AppendTo[nnn,FullSimplify[nnn00]]];
nnn=DeleteDuplicates[nnn]
,
nnn={cor\[Transpose] . nn}];
AppendTo[NNA,nnn]
];AppendTo[NN,NNA]];
tuples={};
For[i1=1,i1<=Length[NN],i1++,
For[i2=1,i2<=Length[NN[[i1]]],i2++,
AppendTo[tuples,Range[Length[NN[[i1,i2]]]]]
]];
tuples=Tuples[tuples];
tmptuples=tuples;
tuples={};
For[i=1,i<=Length[tmptuples],i++,
s3=0;
tuplesA={};
For[i1=1,i1<=Length[NN],i1++,
tuplesB={};
For[i2=1,i2<=Length[NN[[i1]]],i2++,
s3++;
AppendTo[tuplesB,tmptuples[[i,s3]]]];AppendTo[tuplesA,tuplesB]];
AppendTo[tuples,tuplesA]];
RESULT={};
For[i=1,i<=Length[tuples],i++,
tmp1101=hsprelatedV2inlowmsgCor[[imsg,ilow,iset]];
For[i1=1,i1<=Length[hsp],i1++,
For[i2=1,i2<=Length[tmp1101[[i1]]],i2++,
For[i3=1,i3<=Length[tmp1101[[i1,i2,so+1]]],i3++,
tmp1101A=NN[[i1,tmp1101[[i1,i2,so+1,i3,1,1]],tuples[[i,i1,tmp1101[[i1,i2,so+1,i3,1,1]]]]]];
tmp1101B=NN[[i1,tmp1101[[i1,i2,so+1,i3,1,2]],tuples[[i,i1,tmp1101[[i1,i2,so+1,i3,1,2]]]]]];
tmp1101C=tmp1101[[i1,i2,so+1,i3,2]];
If[FullSimplify[tmp1101C\[Transpose].tmp1101A==tmp1101B],Continue[]];
Goto[out1101]
]]];Label[out1101];
If[i1<=Length[hsp],Continue[]];
f1025=0;
bnhsp=bnhsl=bnhspl={};
For[i1=1,i1<=Length[hsp],i1++,
For[i2=1,i2<=Length[hsp[[i1]]],i2++,
nn=NN[[i1,i2,tuples[[i,i1,i2]]]];
For[i3=1,i3<=Length[nn],i3++,
If[IntegerQ[nn[[i3]]],Continue[]];
f1025++;
If[Length[hsp[[i1,i2,6]]]==0,Continue[]];
n3=Position[hsp[[i1,i2,6]]\[Transpose][[1]],{i3}];
If[Length[n3]==0,Continue[]];
AppendTo[bnhsp,{hsp2[[i1,i2,2,3]],hsp2[[i1,i2,2,{2,3}]],{i3},hsp[[i1,i2,6,n3[[1,1]],{2,3}]]}]]
]];
For[i1=1,i1<=Length[hsl],i1++,
For[i2=1,i2<=Length[hsl[[i1]]],i2++,
For[i3=1,i3<=Length[hsl[[i1,i2]]]-1,i3++,
If[Length[hsl[[i1,i2,i3,3]]]==0,Continue[]];
n1=Position[hsp0,hsl[[i1,i2,i3,1,1]]];
s1++;If[Length[n1]==1,s2++];
nn1=NN[[n1[[1,1]],n1[[1,2]],tuples[[i,n1[[1,1]],n1[[1,2]]]]]];
n2=Position[hsp0,hsl[[i1,i2,i3+1,1,1]]];
s1++;If[Length[n2]==1,s2++];
nn2=NN[[n2[[1,1]],n2[[1,2]],tuples[[i,n2[[1,1]],n2[[1,2]]]]]];
If[MemberQ[IntegerQ/@nn1,False],Continue[]];
If[MemberQ[IntegerQ/@nn2,False],Continue[]];
nn1=hsl[[i1,i2,i3,2]]\[Transpose] . nn1;
nn2=hsl[[i1,i2,i3+1,2]]\[Transpose] . nn2;
\[Delta]nn12=nn1-nn2;
n7=n8={};
For[i4=1,i4<=Length[\[Delta]nn12],i4++,
If[\[Delta]nn12[[i4]]>0,AppendTo[n7,i4]];
If[\[Delta]nn12[[i4]]<0,AppendTo[n8,i4]]];
For[i4=1,i4<=Length[n7],i4++,
For[i5=1,i5<=Length[n8],i5++,
n3=Position[hsl[[i1,i2,i3,3]]\[Transpose][[1]],{n7[[i4]],n8[[i5]]}];
If[Length[n3]==1,
AppendTo[bnhsl,{hsl2[[i1,i2,1,2]],{{hsl2[[i1,i2,2,i3,1]],hsl2[[i1,i2,2,i3+1,1]]},{n7[[i4]],n8[[i5]]}},-\[Delta]nn12[[n7[[i4]]]]\[Delta]nn12[[n8[[i5]]]],hsl[[i1,i2,i3,3]][[n3[[1,1]],{2,3}]]}]];
n3=Position[hsl[[i1,i2,i3,3]]\[Transpose][[1]],{n8[[i5]],n7[[i4]]}];
If[Length[n3]==1,
AppendTo[bnhsl,{hsl2[[i1,i2,1,2]],{{hsl2[[i1,i2,2,i3,1]],hsl2[[i1,i2,2,i3+1,1]]},{n8[[i5]],n7[[i4]]}},-\[Delta]nn12[[n7[[i4]]]]\[Delta]nn12[[n8[[i5]]]],hsl[[i1,i2,i3,3]][[n3[[1,1]],{2,3}]]}]
]]];
]]];


For[i1=1,i1<=Length[hspl],i1++,
For[i2=1,i2<=Length[hspl[[i1]]],i2++,
For[i3=1,i3<=Length[hspl[[i1,i2]]],i3++,
For[i33=1,i33<=Length[hspl[[i1,i2]]],i33++,
If[Length[hspl[[i1,i2,i3,3]]]==0,Continue[]];
n1=Position[hsp0,hspl[[i1,i2,i3,1,1]]];
s1++;If[Length[n1]==1,s2++];
nn1=NN[[n1[[1,1]],n1[[1,2]],tuples[[i,n1[[1,1]],n1[[1,2]]]]]];
n2=Position[hsp0,hspl[[i1,i2,i33,1,1]]];
s1++;If[Length[n2]==1,s2++];
nn2=NN[[n2[[1,1]],n2[[1,2]],tuples[[i,n2[[1,1]],n2[[1,2]]]]]];
If[MemberQ[IntegerQ/@nn1,False],Continue[]];
If[MemberQ[IntegerQ/@nn2,False],Continue[]];
nn1=hspl[[i1,i2,i3,2]]\[Transpose] . nn1;
nn2=hspl[[i1,i2,i33,2]]\[Transpose] . nn2;
\[Delta]nn12=nn1-nn2;
n7=n8={};
For[i4=1,i4<=Length[\[Delta]nn12],i4++,
If[\[Delta]nn12[[i4]]>0,AppendTo[n7,i4]];
If[\[Delta]nn12[[i4]]<0,AppendTo[n8,i4]]];
For[i4=1,i4<=Length[n7],i4++,
For[i5=1,i5<=Length[n8],i5++,
n3=Position[hspl[[i1,i2,i3,3]]\[Transpose][[1]],{n7[[i4]],n8[[i5]]}];
If[Length[n3]==1,
AppendTo[bnhspl,{hspl2[[i1,i2,1,2]],{{hspl2[[i1,i2,2,i3,1]],hspl2[[i1,i2,2,i33,1]]},{n7[[i4]],n8[[i5]]}},-\[Delta]nn12[[n7[[i4]]]]\[Delta]nn12[[n8[[i5]]]],hspl[[i1,i2,i3,3]][[n3[[1,1]],{2,3}]]}]];
n3=Position[hspl[[i1,i2,i3,3]]\[Transpose][[1]],{n8[[i5]],n7[[i4]]}];
If[Length[n3]==1,
AppendTo[bnhspl,{hspl2[[i1,i2,1,2]],{{hspl2[[i1,i2,2,i3,1]],hspl2[[i1,i2,2,i33,1]]},{n8[[i5]],n7[[i4]]}},-\[Delta]nn12[[n7[[i4]]]]\[Delta]nn12[[n8[[i5]]]],hspl[[i1,i2,i3,3]][[n3[[1,1]],{2,3}]]}]]
]];
]]]];
If[Length[bnhsp]==0&&Length[bnhsl]==0&&Length[bnhspl]==0&&f1025==0,
AIB=allAIB[[allimsg2[[imsg,ilow]],1,so+1]];
nn={};
For[ii=1,ii<=Length[indhspname[[allimsg2[[imsg,ilow]],1]]],ii++,
n3=Position[hsp00,indhspname[[allimsg2[[imsg,ilow]],1,ii]]];
s1++;If[Length[n3]==1,s2++];
AppendTo[nn,NN[[n3[[1,1]],n3[[1,2]],tuples[[i,n3[[1,1]],n3[[1,2]]]]]]]];
nn=Flatten[nn];
{{dBS,XBS},basis}=AIB;
basis=basis\[Transpose];
q=FullSimplify[Inverse[basis\[Transpose] . basis] . basis\[Transpose] . nn];
s1++;If[DeleteDuplicates[IntegerQ/@FullSimplify[nn]]=={True}&&FullSimplify[basis . q==nn]&&DeleteDuplicates[Table[IntegerQ[FullSimplify[q[[i4]]XBS[[i4]]]],{i4,Length[XBS]}]]=={True},
s2++;
If[DeleteDuplicates[Table[IntegerQ[FullSimplify[q[[i4]]]],{i4,Length[XBS]}]]=={True},case="I",case="II";
SI={};
For[i4=1,i4<=dBS,i4++,
If[XBS[[i4]]>1,AppendTo[SI,Mod[FullSimplify[q[[i4]]XBS[[i4]]],XBS[[i4]]]]]]]];
AppendTo[RESULT,If[case=="I","I",{"II",SI}]],
AppendTo[RESULT,{"III",bnhsp,bnhsl,bnhspl}]
]];
s1++;If[Length[RESULT]>0,s2++];
{{s1,s2},RESULT}]


SIValues[imsg_,so_,SI_]:=Module[{XBS,i,j,re},XBS=allAIB[[imsg,1,so+1,1,2]];
re={};
For[i=1,i<=Length[XBS],i++,
If[XBS[[i]]>1,
AppendTo[re,XBS[[i]]]
]];
{SI,Subscript["Z",ToString[re]]}
]
ShowBC[result_]:=Table[{result[[i1,1]],result[[i1,2]],result[[i1,3]],{result[[i1,4,1]],MatrixForm[result[[i1,4,2]]]}},{i1,Length[result]}]


CalTopoEvol[msg_,symdata_,so_]:=Module[{topoevol,sub,patten,path,path1,patten2,allimsg3,ipath,imsg},
imsg=Position[allmsg3,msg][[1,1]];
topoevol=Table[CalTopoSubGrp[imsg,SBPattern[[imsg,ii]],symdata,so],{ii,Length[SBPattern[[imsg]]]}];
sub=Position[allmpg,mpgmsg[[imsg]]][[1,1]];
patten=SBPattern[[imsg]];
path=resultfrom1toi[[sub]];
sub=submpg[[sub,2]];
path1=Table[{patten[[Position[sub,path[[i,j]]][[1,1]]]],patten[[Position[sub,path[[i,j+1]]][[1,1]]]]},{i,Length[path]},{j,Length[path[[i]]]-1}];ipath=Input["Specify symmetry-breaking paths ("<>ToString[Length[path]]<>" paths in total): e.g. {1,2,3}"];patten2=DeleteDuplicates[Flatten[path1[[ipath]]]];
allimsg3=Table[lowmsg[[imsg,Position[patten,patten2[[iii]]][[1,1]],8]],{iii,Length[patten2]}];
RelationGraph[Length[Position[path1[[ipath]],{#2,#1}]]>0&,patten2,patten2,GraphLayout->"LayeredDigraphEmbedding",VertexLabels->Table[patten2[[iii]]->Placed[PopupWindow[ToString[patten2[[iii]]],Grid[Flatten[{{{"{Symmetry-breaking pattern,t-subgroup}",{patten2[[iii]],allimsg3[[iii]]},SpanFromLeft}},Table[If[Length[topoevol[[Position[patten,patten2[[iii]]][[1,1]],2]][[iiii]]]==0,{topoevol[[Position[patten,patten2[[iii]]][[1,1]],2]][[iiii]],"Trivial by symmetry-indicators",SpanFromLeft},

If[topoevol[[Position[patten,patten2[[iii]]][[1,1]],2]][[iiii,1]]=="II",
{topoevol[[Position[patten,patten2[[iii]]][[1,1]],2]][[iiii,1]],SIValues[Position[allmsg3,allimsg3[[iii]]][[1,1]],1,topoevol[[Position[patten,patten2[[iii]]][[1,1]],2]][[iiii,2]]],SpanFromLeft},
{topoevol[[Position[patten,patten2[[iii]]][[1,1]],2]][[iiii,1]],
If[Length[topoevol[[Position[patten,patten2[[iii]]][[1,1]],2]][[iiii,2]]]==0,"",Grid[ShowBC[topoevol[[Position[patten,patten2[[iii]]][[1,1]],2]][[iiii,2]]],Frame->All,FrameStyle->Orange]],
If[Length[topoevol[[Position[patten,patten2[[iii]]][[1,1]],2]][[iiii,3]]]==0,"",Grid[ShowBC[topoevol[[Position[patten,patten2[[iii]]][[1,1]],2]][[iiii,3]]],Frame->All,FrameStyle->Orange]],
If[Length[topoevol[[Position[patten,patten2[[iii]]][[1,1]],2]][[iiii,4]]]==0,"",PopupWindow["Click to show details",Grid[ShowBC[topoevol[[Position[patten,patten2[[iii]]][[1,1]],2]][[iiii,4]]],Frame->All,FrameStyle->Orange],
WindowElements->{"VerticalScrollBar", "HorizontalScrollBar"}
]]
}
]
],{iiii,Length[topoevol[[Position[patten,patten2[[iii]]][[1,1]],2]]]}]},1],Frame->All],WindowSize->Full,WindowTitle->"Topological states based on the t-subgroup, "<>ToString[{patten2[[iii]],allimsg3[[iii]]}],
WindowElements->{"VerticalScrollBar", "HorizontalScrollBar"}
],Center],{iii,Length[patten2]}],
VertexStyle->Directive[EdgeForm[Red],White],EdgeStyle->Blue,VertexSize->Large
]]


End[];


EndPackage[];
