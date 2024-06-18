(* ::Package:: *)

(* ::Title:: *)
(*Band Utilities*)


(* ::Author:: *)
(*Anirudh Chandrasekaran*)


(* ::Affiliation:: *)
(*8 February 2024*)


(* ::Text:: *)
(*We include routines that are needed for plotting band-structure, writing Hamiltonians to C++ compatible files, and fitting Hamiltonians to bands along k-paths and generating Fermi surfaces by brute force. The Fermi-surface functionality is yet to be refined.*)


(* ::Section::Closed:: *)
(*Projecting to orbitals, i.e finding the orbital/orbital group with maximum contribution*)


FindDominantOrbital[EigVec_,OrbitalGroups_,OrbitalTrfmn_,Resln_:10^-6]:=
Module[{EigProjn},
EigProjn=First@Last@Reap[Do[
Sow[Total[Abs[(OrbitalTrfmn . EigVec)[[Grp]]]^2]]
,{Grp,OrbitalGroups}]];

RandomChoice[Flatten[Position[Chop[EigProjn-Max[EigProjn],Resln],0]]]
]


(* ::Section::Closed:: *)
(*Checking the type and validity of a piecewise linear path in the Brillouin zone specified by the end points of segments*)


CheckPathType[HighSymmPtsList_]:=
Module[
{SpaceDim},
If[!ListQ[HighSymmPtsList]||!AllTrue[HighSymmPtsList,ListQ[#]&]||!AllTrue[HighSymmPtsList,ListQ[#[[1]]]&],
Return[{"Invalid",0}]
,

SpaceDim=Length[HighSymmPtsList[[1,1]]];
If[AllTrue[HighSymmPtsList,Length[#[[1]]]==SpaceDim&] && AllTrue[Flatten[HighSymmPtsList[[All,1]]]//N,NumberQ[#]&] && AllTrue[Flatten[HighSymmPtsList[[All,2]]],StringQ[#]&],
Return[{"Single",SpaceDim}]
,

SpaceDim=Length[HighSymmPtsList[[1,1,1]]];

Do[
If[SpaceDim==0||!AllTrue[HighSymmPts,Length[#[[1]]]==SpaceDim&] || !AllTrue[Flatten[HighSymmPts[[All,1]]]//N,NumberQ[#]&] || !AllTrue[Flatten[HighSymmPts[[All,2]]],StringQ[#]&],
Return[{"Invalid",0}]
]
,{HighSymmPts,HighSymmPtsList}];

Return[{"Multi",SpaceDim}]
];
];

]


(* ::Section::Closed:: *)
(*Function to generate the band-structure data*)


(* ::Text:: *)
(*The following function generates the band-structure of a Hamiltonian along a k-path in the Brillouin Zone (BZ). Two possible choices of k-path are allowed: *)
(*1) a continuous, piece-wise linear path connecting some chosen points in the BZ (possibly, but not necessarily high symmetry points) and *)
(*2) a discontinuous path composed of a collection of piece-wise linear continuous paths passing through some specified points. *)
(*To give an example, in the former case, we could feed as input the following list: {{{0, 0}, "\[CapitalGamma]"}, {{\[Pi], \[Pi]}, "M"}, {{\[Pi], 0}, "X"}, {{0, 0}, "\[CapitalGamma]"}}, which is an ordered collection of the coordinates of the points and their symbols. *)
(*In the latter case, we could give two or more continuous, piece-wise linear sub-paths, for example: {{{{0, 0}, "\[CapitalGamma]"}, {{\[Pi], \[Pi]}, "M"}, {{\[Pi], 0}, "X"}, {{0, 0}, "\[CapitalGamma]"}}, {{{-\[Pi], 0}, "Subscript[X, 1]"}, {{0, \[Pi]}, "Subscript[X, 2]"}}}.*)
(**)
(*The final output is a collection of points, each with a one dimensional k-coordinate (length) and the set of all the band eigenvalues at that point. We can plot this directly since the k-coordinate is adjusted to enable plotting. Alternatively we can also output this to a file for plotting with another utility like Matplotlib or GNU Plot. Other useful data like the locations of the appropriate labels (of high-symmetry points) and the total number of points are also given as a part of the output.*)
(**)


Options[GenerateBandStructure]={"TotalPoints"->200,"Resolution"->10^-6,"OrbitalProjection"->False,"OrbitalGrouping"->-1,"OrbitalRedefinition"->-1}

GenerateBandStructure[Hamltn_,HighSymmPtsList_,OptionsPattern[]]:=
Module[{DirVecs,SpaceDim,TotalPts,TotalPtsIn,ReslnIn,Resln,TotalLen,PathType,NoPtsSeg,LenSeg,WithProjections,TotalPtsSubPath,TotalPtsFinal,OrbitalGroups,OrbitalTrfmn,HighSymmPts,BandStructure,BandStructureList,BandFull,pVec,CurPt,CurLen,HamltnDim,LabelLocs,LabelLocsLength,LabelLocsPos,LenSubPath,CurSubPath,AllFine,EigsSys,BandsProj},

TotalPtsIn=OptionValue["TotalPoints"];
ReslnIn=OptionValue["Resolution"];


(*We begin with some sanity checks on the last two inputs, namely the number of plot points and the eigenvalue resolution.*)
If[IntegerQ[TotalPtsIn]&&TotalPtsIn>0,
TotalPts=TotalPtsIn;
,
Print[Style["Invalid number of plot points provided! Reverting to the default value of 200.",Red]];
TotalPts=200;
];

If[NumberQ[ReslnIn]&&ReslnIn\[Element]Reals&&ReslnIn>0,
Resln=ReslnIn;
,
Print[Style["Invalid eigenvalue-resolution provided! Reverting to the default value of \!\(\*SuperscriptBox[\(10\), \(-6\)]\).",Red]];
Resln=10^-6;
];

{PathType,SpaceDim}=CheckPathType[HighSymmPtsList];

(*To compute the dimension of the Hamiltonian matrix, we feed use the H-matrix evaluated at the origin:*)
HamltnDim=Length[Hamltn[ConstantArray[0,SpaceDim]]];

(*Print[HamltnDim];*)

(*Some sanity checks on the Hamiltonian.*)
If[!(HamltnDim>=1 && AllTrue[AllTrue[NumericQ]][Hamltn[ConstantArray[0,SpaceDim]]//N]),
Print[Style["Invalid Hamiltonian provided. The Hamiltonian can not contain symbolic terms when evaluated at a numerical \!\(\*
StyleBox[\"k\",\nFontSlant->\"Italic\"]\)-vector!",Red]];
Return[Null]
];

(*We now work out the orbital projection scheme based on the input options*)
OrbitalGroups={};

(*Print[OptionValue["OrbitalProjection"]];*)

If[OptionValue["OrbitalProjection"]==True,

If[!ListQ[OptionValue["OrbitalRedefinition"]]||!AllTrue[Flatten[OptionValue["OrbitalRedefinition"]],NumberQ[#]&]||Dimensions[OptionValue["OrbitalRedefinition"]]!={HamltnDim,HamltnDim},
OrbitalTrfmn=IdentityMatrix[HamltnDim];
,
OrbitalTrfmn=OptionValue["OrbitalRedefinition"];
];

If[ListQ[OptionValue["OrbitalGrouping"]] && AllTrue[Flatten[OptionValue["OrbitalGrouping"]],IntegerQ[#] && 1<=#<=HamltnDim&],

If[Length[Complement[Table[i,{i,HamltnDim}],Flatten[OptionValue["OrbitalGrouping"]]]]==0,
OrbitalGroups=First@Last@Reap[Do[
If[ListQ[elm],
Sow[elm]
,
Sow[{elm}]
]
,{elm,OptionValue["OrbitalGrouping"]}]];
,
OrbitalGroups=Append[OptionValue["OrbitalGrouping"],Complement[Table[i,{i,HamltnDim}],Flatten[OptionValue["OrbitalGrouping"]]]];
]
,
OrbitalGroups=Table[{i},{i,HamltnDim}];
];

WithProjections=True;
(*Print[OrbitalGroups];*)
,
WithProjections=False;
];

(*We now fork out depending on whether the path is continuous or contains a number of continuous segments with possible discontinuities in between.*)
If[PathType=="Single",(*This means that the k-path provided contains a single, piece-wise linear and continuous segment.*)
(*Length[Dimensions[HighSymmPtsList]]==2*)
(*Print["Single Path!"];*)

SpaceDim=Length[HighSymmPtsList[[1,1]]];
HighSymmPts=HighSymmPtsList;

(*DirVecs is a list that contains the direction vectors connecting every two consecutive high-symmetry points of each linear segment in the k-path.*)
(*The number of such vectors is one less than the number of high-symmetry points on the k-path.*)
DirVecs=First@Last@Reap[Do[
Sow[HighSymmPts[[i+1,1]]-HighSymmPts[[i,1]]]
,{i,1,Length[HighSymmPts]-1}
]
];


(*We will now calculate the total k-length of the k-path below, after some preliminary calculations.*)
TotalLen=0; 

(*A list that will contain the k-length of each continuous linear segment of the k-path:*)
LenSeg=First@Last@Reap[Do[
Sow[Norm[dir]]
,{dir,DirVecs}
]
];

Do[TotalLen+=elem,{elem ,LenSeg}];


(*The following list (NoPtsSeg) will contain the number of plot points in each segment of the k-path, such that the k-distance between adjacent points is constant throughout.*)
(*The loop to compute this list NoPtsSeg, is to be executed AFTER the loop to compute TotalLen since we need to first compute TotalLen, which is the used to compute the number of points in each segment, based on the ratio between the length of the segment and the total length*)
NoPtsSeg=First@Last@Reap[Do[
Sow[Round[1.*TotalPts*Norm[dir]/TotalLen]]
,{dir,DirVecs}
]
];

(*We need to recalculate total no of points since it may have changed by an integer close to 1*)
TotalPtsFinal=0; 
Do[
TotalPtsFinal+=elem
,{elem ,NoPtsSeg}
];

(*Finally we go on to compute the bandstructure.*)
BandStructure=ConstantArray[0,{TotalPtsFinal,1}];
CurLen=0;
CurPt=0; (*The index number of the first k-point in the current linear segment of the k-path, for which the bandstructure is being calculated.*)
pVec=ConstantArray[0,SpaceDim];

If[WithProjections==True,
Do[(*iterating through the segments*)
Do[
(*The first element is the point's location measured in terms of k-length, used as a coordinate.*)
BandStructure[[CurPt+j,1]]=CurLen+1.*(j-1.)*LenSeg[[i]]/NoPtsSeg[[i]]; 

pVec=HighSymmPts[[i,1]]+1.*(j-1.)*DirVecs[[i]]/NoPtsSeg[[i]];

EigsSys=Eigensystem[0.5(N[Hamltn[pVec]]+ConjugateTranspose[N[Hamltn[pVec]]])]//N;
(*To ensure that the Hamiltonian is Hermitian, we take an average of H and H^\[Dagger]*)
BandsProj=Chop[Table[{EigsSys[[1,l]],FindDominantOrbital[EigsSys[[2,l]],OrbitalGroups,OrbitalTrfmn,Resln]},{l,HamltnDim}],Resln];
(*Ordering[Abs[EigsSys[[2,l]]],-1][[1]]*)

BandStructure[[CurPt+j]]=Flatten[Join[BandStructure[[CurPt+j]],SortBy[BandsProj,#[[1]]&]]];

,{j,1,NoPtsSeg[[i]]}];

CurPt=CurPt+NoPtsSeg[[i]];
CurLen=CurLen+LenSeg[[i]];

,{i,1,Length[HighSymmPts]-1}];

,
(*No orbital projections*)
Do[(*iterating through the segments*)
Do[
(*The first element is the point's location measured in terms of k-length, used as a coordinate.*)
BandStructure[[CurPt+j,1]]=CurLen+1.*(j-1.)*LenSeg[[i]]/NoPtsSeg[[i]]; 

pVec=HighSymmPts[[i,1]]+1.*(j-1.)*DirVecs[[i]]/NoPtsSeg[[i]];
BandStructure[[CurPt+j]]=Join[BandStructure[[CurPt+j]],Sort[Chop[Eigenvalues[0.5(N[Hamltn[pVec]]+ConjugateTranspose[N[Hamltn[pVec]]])]//N,Resln]]]
(*To ensure that the Hamiltonian is Hermitian, we take an average of H and H^\[Dagger]*)

,{j,1,NoPtsSeg[[i]]}];

CurPt=CurPt+NoPtsSeg[[i]];
CurLen=CurLen+LenSeg[[i]];

,{i,1,Length[HighSymmPts]-1}];
];

(*Finally we find the x-coordinate where the symbols/labels of the high symmetry points are to be placed in the final plot:*)
LabelLocs=ConstantArray[0,{Length[HighSymmPts],2}];
LabelLocs[[1,2]]=HighSymmPts[[1,2]];

Do[
LabelLocs[[i,1]]=LabelLocs[[i-1,1]]+LenSeg[[i-1]]; (*or = NoPtsSeg[[i-1]] if coordinate is to be given in the point-number system rather than in terms of path-length system.*)
LabelLocs[[i,2]]=HighSymmPts[[i,2]];
,{i,2,Length[HighSymmPts]}];

(*The list we shall return will contain the bandstructure data-table, location of labels, dimension of Hamiltonian, dimension of space and total number of plot points. This is all the information that is needed for plotting the bandstructure.*)
BandFull={BandStructure,LabelLocs,HamltnDim,SpaceDim,TotalPtsFinal,TotalLen,OptionValue["OrbitalProjection"],OrbitalGroups,PathType};
BandFull

,

If[PathType=="Multi", (*This means that the path contains at least two (or more) piece-wise linear & continuous segments with a possible discontinuity in between these.*)
(*Length[Dimensions[HighSymmPtsList]]==1*)
(*Print["Multi Path!"];*)

LabelLocsLength=1; (*This is the length (or number of elements) of the list LabelLocs that contains the locations and labels of high-symmetry points. At discontinuities we have one label representing two high-symmetry points for example \[CapitalGamma]/Z. After the loop below we would be left with a value one less than the actual length. So we start here with 1 to compensate*)

TotalLen=0; (*Total k-length of the plot*)

TotalPtsFinal=0;

LenSubPath=First@Last@Reap[
(*LenSubPath[[i]] tells us the combined k-length of all the preceeding subpaths, i.e up to and including the (i-1)^th subpath. (Recall that we have a number of piecewise linear subpaths constituting the full k-path). The reason we compute this is that it is simply used as an offset later, while generating the actual dataset. We shall be generating the data separately for each of the subaths and this list contains the starting k values of the subpaths.*)
Do[
(*We iterate through each k-subpath in HighSymmPtsList, and sum over norms of the connecting vectors, to compute TotalLen. The number of such connecting vectors helps determine the number of high-symmetry point labels.*)
Sow[TotalLen];
Do[
TotalLen+=Norm[HighSymmPts[[i+1,1]]-HighSymmPts[[i,1]]];
LabelLocsLength+=1;
,{i,1,Length[HighSymmPts]-1}
]
,{HighSymmPts,HighSymmPtsList}] (*Each element of HighSymmPtsList is a continuous subpath of the contour*)
];

CurSubPath=1; (*An index to keep track of the current k-subpath for which the band structure is being generated*)

BandStructureList=First@Last@Reap[
Do[ (*Here again we iterate through each subpath defined by the element HighSymmPts contained in the list called HighSymmPtsList*)

(*We first compute the list containing the connecting vectors in the piecewise linear k-subpath*)
DirVecs=First@Last@Reap[Do[
Sow[HighSymmPts[[i+1,1]]-HighSymmPts[[i,1]]]
,{i,1,Length[HighSymmPts]-1}
]
];

(*k-length of each linear segment in the k-subpath*)
LenSeg=First@Last@Reap[Do[
Sow[Norm[dir]]
,{dir,DirVecs}
]
];

(*number of plot points in each linear segment in the k-subpath*)
NoPtsSeg=First@Last@Reap[Do[
Sow[Round[1.*TotalPts*Norm[dir]/TotalLen]]
,{dir,DirVecs}
]
];


(*We need to calculate total no of points in each subpath. Note that this is different from the variable TotalPts which is the full total number of points in the plot*)
TotalPtsSubPath=0; 
Do[
TotalPtsSubPath+=elem,{elem ,NoPtsSeg}
];

TotalPtsFinal+=TotalPtsSubPath;

BandStructure=ConstantArray[0,{TotalPtsSubPath,1}];

CurLen=LenSubPath[[CurSubPath]]; (*This ensures that the different segments of the band-structure corresponding to the different subpaths have the right k-length offset for the first point*)
CurPt=0; (*The index number of the k-point of the plot, for which the bandstructure is currently being calculated*)
pVec=ConstantArray[0,SpaceDim];

If[WithProjections==True,
Do[
Do[
BandStructure[[CurPt+j,1]]=CurLen+1.*(j-1.)*LenSeg[[i]]/NoPtsSeg[[i]];
pVec=HighSymmPts[[i,1]]+1.*(j-1.)*DirVecs[[i]]/NoPtsSeg[[i]];

EigsSys=Eigensystem[0.5(N[Hamltn[pVec]]+ConjugateTranspose[N[Hamltn[pVec]]])]//N;
BandsProj=Chop[Table[{EigsSys[[1,l]],FindDominantOrbital[EigsSys[[2,l]],OrbitalGroups,OrbitalTrfmn,Resln]},{l,HamltnDim}],Resln];
(*Ordering[Abs[EigsSys[[2,l]]],-1][[1]]*)

BandStructure[[CurPt+j]]=Flatten[Join[BandStructure[[CurPt+j]],SortBy[BandsProj,#[[1]]&]]];
,{j,1,NoPtsSeg[[i]]}];
CurPt=CurPt+NoPtsSeg[[i]];
CurLen=CurLen+LenSeg[[i]];
,{i,1,Length[HighSymmPts]-1}];

,
(*No orbital projection*)
Do[
Do[
BandStructure[[CurPt+j,1]]=CurLen+1.*(j-1.)*LenSeg[[i]]/NoPtsSeg[[i]];
pVec=HighSymmPts[[i,1]]+1.*(j-1.)*DirVecs[[i]]/NoPtsSeg[[i]];
BandStructure[[CurPt+j]]=Join[BandStructure[[CurPt+j]],Sort[Chop[Eigenvalues[0.5(N[Hamltn[pVec]]+ConjugateTranspose[N[Hamltn[pVec]]])]//N,Resln]]];
,{j,1,NoPtsSeg[[i]]}];
CurPt=CurPt+NoPtsSeg[[i]];
CurLen=CurLen+LenSeg[[i]];
,{i,1,Length[HighSymmPts]-1}];

];

CurSubPath+=1;

Sow[BandStructure]

,{HighSymmPts,HighSymmPtsList}]
]; 

LabelLocs=ConstantArray[0,{LabelLocsLength,2}];
LabelLocs[[1,2]]=HighSymmPtsList[[1]][[1,2]];
LabelLocsPos=1;
CurSubPath=1;
Do[
If[LabelLocsPos!=1,(*The first symbol of any subpath (other than the first) has to be joined with the last symbol of the preceeding subpath along with a "/". For example, \[CapitalGamma]/M would be the symbol at the interface of two disconnected subpaths, the first ending at \[CapitalGamma] and the second beginning at M*)
LabelLocs[[LabelLocsPos,2]]=LabelLocs[[LabelLocsPos,2]]<>"/"<>HighSymmPts[[1,2]]
];

Do[
LabelLocsPos+=1;
LabelLocs[[LabelLocsPos,1]]=LabelLocs[[LabelLocsPos-1,1]]+Norm[HighSymmPts[[i,1]]-HighSymmPts[[i-1,1]]];

LabelLocs[[LabelLocsPos,2]]=HighSymmPts[[i,2]];
,{i,2,Length[HighSymmPts]}];

CurSubPath+=1;

,{HighSymmPts,HighSymmPtsList}
];

BandFull={BandStructureList,LabelLocs,HamltnDim,SpaceDim,TotalPtsFinal,TotalLen,OptionValue["OrbitalProjection"],OrbitalGroups,PathType};
BandFull

,

Print[Style["Incorrect format used for the list containing high-symmetry points! Execute GenerateBandStructure::usage to see the correct usage.",Red]]
]
]
]



(* ::Section::Closed:: *)
(*Plotting the generated band-structure*)


Options[PlotBandStructure]={"AspectRatio"->0.75, "xLabel"->{(*Font size*)12,(*Font colour*)Black,(*Font family*)"Arial"}, "yLabel"->{(*Label*)"E",(*Font Size*)12,(*Font Color*)Black,(*Font Family*)"Arial"}, "yTicks"->{10,Black,"Helvetica"}, "LineThickness"->0.005, "LineColorScheme"->{Blue}, "DividingLines"->{(*Dashing[..]*)0.02,(*Color*)Gray,(*Thickness*)0.005}, "yLabelRotate"->0, "PlotKeyLegend"->{None,{0,0}}}

PlotBandStructure[BandStructure_,yRange_,OptionsPattern[]]:=
Module[{DivLines,DivLinesGraphics,EpilogPlot,NullTable,xLabelFont,yLabel,PlotAspectRatio,yTicsFont,LineThickness,LineColorScheme,DivLineScheme,yLabelRotate,PlotKeyLegend,BandPlot,OrbitalColorList,BandData,BandList,CurSubPath,LineStyle},

PlotAspectRatio=OptionValue["AspectRatio"];
xLabelFont=OptionValue["xLabel"];
yLabel=OptionValue["yLabel"];
yTicsFont=OptionValue["yTicks"];
LineThickness=OptionValue["LineThickness"];
(*LineColorScheme=OptionValue["LineColorScheme"];*)
DivLineScheme=OptionValue["DividingLines"];
yLabelRotate=OptionValue["yLabelRotate"];
PlotKeyLegend=OptionValue["PlotKeyLegend"];

(*DivLines is a list containing the plot directives for the dashed, dividing vertical lines that separate the different segments in the sphagetti bandstructure. These are located at the plot points corresponding to the high-symmetry points making up the piecewise linear k-path*)
DivLines=ConstantArray[0,Length[BandStructure[[2]]]-2];(*BandStructure[[2]] is the list LabelLocs containing the x-coordinates of the labels of the high symmetry points*)
Do[
DivLines[[i]]={Dashing[DivLineScheme[[1]]],DivLineScheme[[2]],Thickness[DivLineScheme[[3]]],Line[{{BandStructure[[2]][[i+1,1]],yRange[[1]]},{BandStructure[[2]][[i+1,1]],yRange[[2]]}}]}(*BandStructure[[2]] is LabelLocs*)
,{i,1,Length[DivLines]}];

If[BandStructure[[7]]==False,
If[ColorQ[OptionValue["LineColorScheme"]],
LineColorScheme={OptionValue["LineColorScheme"]};
,
If[ListQ[OptionValue["LineColorScheme"]]&&AllTrue[OptionValue["LineColorScheme"],ColorQ[#]&],
LineColorScheme=OptionValue["LineColorScheme"];
,
LineColorScheme={Blue};
]
];

LineStyle=Table[{LineColorScheme[[Mod[i-1,Length[LineColorScheme]]+1]],Thickness[LineThickness]},{i,1,BandStructure[[3]]}]; (*BandStructure[[3]] is the dimension of the Hamiltonian.*)
,
If[!ListQ[OptionValue["LineColorScheme"]]||!AllTrue[OptionValue["LineColorScheme"],ColorQ[#]&]||Length[OptionValue["LineColorScheme"]]!=Length[BandStructure[[8]]],
LineColorScheme=RandomColor[Length[BandStructure[[8]]]];
,
LineColorScheme=OptionValue["LineColorScheme"];
];

LineStyle=Table[{LineColorScheme[[Mod[i-1,Length[LineColorScheme]]+1]],Thickness[LineThickness]},{i,1,Length[BandStructure[[8]]]}]; (*Length[BandStructure[[8]]] is the number of orbital groups.*)
NullTable=Table[Null,{j,Length[LineStyle[[All,1]]]}];

If[!ListQ[OptionValue["PlotKeyLegend"][[1]]]||!AllTrue[OptionValue["PlotKeyLegend"][[1]],StringQ[#]&]||Length[OptionValue["PlotKeyLegend"][[1]]]!=Length[BandStructure[[8]]],
PlotKeyLegend[[1]]=StringRiffle[#,","]&/@BandStructure[[8]];
]
];

(*We use the Graphics command to display these lines together*)
DivLinesGraphics=Graphics[DivLines];


(*This list contains the column index of the bands as stored in the table. This needed below for plotting.*)
BandList=Table[i,{i,2,BandStructure[[3]]+1}];

(*We now display them together along with the user-supplied styling*)


If[BandStructure[[9]]=="Single",

If[BandStructure[[7]]==False,(*That is, no orbital projections*)
BandPlot=ListLinePlot[BandStructure[[1]][[All,{1,#}]]&/@BandList,PlotRange->{{0,BandStructure[[6]]},yRange},AspectRatio->PlotAspectRatio,PlotStyle->LineStyle,Frame->True,FrameTicks->{{Automatic,None},{BandStructure[[2]](*BandStructure[[2]] is LabelLocs*),None}},FrameTicksStyle->{{Directive[FontSize->yTicsFont[[1]],FontColor->yTicsFont[[2]],FontFamily->yTicsFont[[3]]],Automatic},{Directive[FontSize->xLabelFont[[1]],FontColor->xLabelFont[[2]],FontFamily->xLabelFont[[3]]],Automatic}},PlotRangePadding->0,FrameLabel->{{Rotate[Style[yLabel[[1]],FontSize->yLabel[[2]],FontColor->yLabel[[3]],FontFamily->yLabel[[4]]],yLabelRotate],None},{None,None}},
PlotLegends->Placed[{PlotKeyLegend[[1]]},{Scaled[PlotKeyLegend[[2]]],{0,0.5}(*This marks the left edge center of the label*)}]];

Show[BandPlot,DivLinesGraphics]
,

OrbitalColorList=LineStyle[[All,1]][[#]]&/@Flatten/@(BandStructure[[1]][[All,{#}]]&/@Table[2i+1,{i,1,BandStructure[[3]]}]);

EpilogPlot=Plot[NullTable,{x,0,BandStructure[[6]]},PlotStyle->LineStyle[[All,1]],PlotRange->{{0,BandStructure[[6]]},yRange}
,
PlotLegends->Placed[PlotKeyLegend[[1]],{Scaled[PlotKeyLegend[[2]]],{0,0.5}(*This marks the left edge center of the label*)}]
];

BandPlot=
Table[ListLinePlot[BandStructure[[1]][[All,{1,2i}]],Mesh->Length[BandStructure[[1]][[All,{1,2i}]]],MeshStyle->Transparent,MeshShading->OrbitalColorList[[i]],PlotRange->{{0,BandStructure[[6]]},yRange}
,PlotStyle->Thickness[LineThickness],AspectRatio->PlotAspectRatio,Frame->True,FrameTicks->{{Automatic,None},{BandStructure[[2]](*BandStructure[[2]] is LabelLocs*),None}},FrameTicksStyle->{{Directive[FontSize->yTicsFont[[1]],FontColor->yTicsFont[[2]],FontFamily->yTicsFont[[3]]],Automatic},{Directive[FontSize->xLabelFont[[1]],FontColor->xLabelFont[[2]],FontFamily->xLabelFont[[3]]],Automatic}},PlotRangePadding->0,FrameLabel->{{Rotate[Style[yLabel[[1]],FontSize->yLabel[[2]],FontColor->yLabel[[3]],FontFamily->yLabel[[4]]],yLabelRotate],None},{None,None}}
],{i,BandStructure[[3]]}];

Show[BandPlot,EpilogPlot,DivLinesGraphics]
]

,

If[BandStructure[[9]]=="Multi",
CurSubPath=1;
BandPlot=
First@Last@Reap[
Do[
If[CurSubPath==1,(*The first subpath is singled out since it is used to define all the plot features and parameters. The rest of the segments are plotted bare bones and combined with the first segment at the end. This avoids unnecessary complications and conflicts*)

If[BandStructure[[7]]==False,
Sow[ListLinePlot[BandData[[All,{1,#}]]&/@BandList,PlotRange->{{0,BandStructure[[6]]},yRange},AspectRatio->PlotAspectRatio,PlotStyle->LineStyle,Frame->True,FrameTicks->{{Automatic,None},{BandStructure[[2]](*BandStructure[[2]] is LabelLocs*),None}},FrameTicksStyle->{{Directive[FontSize->yTicsFont[[1]],FontColor->yTicsFont[[2]],FontFamily->yTicsFont[[3]]],Automatic},{Directive[FontSize->xLabelFont[[1]],FontColor->xLabelFont[[2]],FontFamily->xLabelFont[[3]]],Automatic}},PlotRangePadding->0,FrameLabel->{{Rotate[Style[yLabel[[1]],FontSize->yLabel[[2]],FontColor->yLabel[[3]],FontFamily->yLabel[[4]]],yLabelRotate],None},{None,None}},PlotLegends->Placed[{PlotKeyLegend[[1]]},{Scaled[PlotKeyLegend[[2]]],{0,0.5}(*This marks the left edge center of the label*)}]]]
,

OrbitalColorList=LineStyle[[All,1]][[#]]&/@Flatten/@(BandData[[All,{#}]]&/@Table[2i+1,{i,1,BandStructure[[3]]}]);

Sow[Table[ListLinePlot[BandData[[All,{1,2i}]],Mesh->Length[BandData[[All,{1,2i}]]],MeshStyle->Transparent,MeshShading->OrbitalColorList[[i]],PlotRange->{{0,BandStructure[[6]]},yRange}
,PlotStyle->Thickness[LineThickness],AspectRatio->PlotAspectRatio,Frame->True,FrameTicks->{{Automatic,None},{BandStructure[[2]](*BandStructure[[2]] is LabelLocs*),None}},FrameTicksStyle->{{Directive[FontSize->yTicsFont[[1]],FontColor->yTicsFont[[2]],FontFamily->yTicsFont[[3]]],Automatic},{Directive[FontSize->xLabelFont[[1]],FontColor->xLabelFont[[2]],FontFamily->xLabelFont[[3]]],Automatic}},PlotRangePadding->0,FrameLabel->{{Rotate[Style[yLabel[[1]],FontSize->yLabel[[2]],FontColor->yLabel[[3]],FontFamily->yLabel[[4]]],yLabelRotate],None},{None,None}}
],{i,BandStructure[[3]]}]];
]
,

If[BandStructure[[7]]==False,
Sow[ListLinePlot[BandData[[All,{1,#}]]&/@BandList,PlotRange->yRange,AspectRatio->PlotAspectRatio,PlotStyle->LineStyle,Frame->None,Ticks->None,Axes->False]]
,

OrbitalColorList=LineStyle[[All,1]][[#]]&/@Flatten/@(BandData[[All,{#}]]&/@Table[2i+1,{i,1,BandStructure[[3]]}]);

Sow[Table[ListLinePlot[BandData[[All,{1,2i}]],Mesh->Length[BandData[[All,{1,2i}]]],MeshStyle->Transparent,MeshShading->OrbitalColorList[[i]],PlotRange->yRange
,PlotStyle->Thickness[LineThickness],AspectRatio->PlotAspectRatio,Frame->None,Ticks->None,Axes->False
],{i,BandStructure[[3]]}]];
]
];

CurSubPath+=1;
,{BandData,BandStructure[[1]]}];

Sow[DivLinesGraphics];
];

If[BandStructure[[7]]==False,
Show[BandPlot]
,

EpilogPlot=Plot[NullTable,{x,0,BandStructure[[6]]},PlotStyle->LineStyle[[All,1]],PlotRange->{{0,BandStructure[[6]]},yRange}
,
PlotLegends->Placed[PlotKeyLegend[[1]],{Scaled[PlotKeyLegend[[2]]],{0,0.5}(*This marks the left edge center of the label*)}]
];

Show[BandPlot,EpilogPlot]
]
]
]
]


(* ::Section::Closed:: *)
(*Diagnosing higher order singularities*)


DiagnoseHOS[Poly_,MaxDeg_:-1,Resln_:10^-6]:=Module[
{PolyDeg,SpaceDim,KVec,CoeffRow,Corank,HessMat,HessEigs,JacobMat,PolyDerivs,Monoms,AllMonoms,Partns,PostvDistPartns,PartnsWithZeros,BasisPolys,CoeffMatrix,ZeroList,OneList,BMat,SolMat,Solved,Determinacy,Codimension},

KVec=Variables[Poly]; (*We detect the symbolic variables used in the polynomial.*)
(*Print[KVec];*)

SpaceDim=Length[KVec];

(*To compute the corank, we need to compute the Hessian matrix first:*)
HessMat=First@Last@Reap[
Do[
Sow[
First@Last@Reap[
Do[
Sow[ReplaceAll[D[Poly,K1,K2],Table[Ki->0,{Ki,KVec}]]];
,{K2,KVec}]
]
]
,{K1,KVec}]
];

HessEigs=Chop[Eigenvalues[HessMat],Resln]//Rationalize;

(*The corank is the number of zero eigenvalues of the Hessian:*)
Corank=Count[HessEigs,0];

(*The Jacobian is necessary to determine if we have a critical point in the first place:*)
JacobMat=Chop[First@Last@Reap[
Do[
Sow[ReplaceAll[D[Poly,K1],Table[Ki->0,{Ki,KVec}]]]
,{K1,KVec}]
],Resln]//Rationalize;

(*Print[JacobMat];*)
Solved=False;

If[Count[JacobMat,0]==SpaceDim,(*If all the elements of the Jacobian are zero, we have a critical point.*)
If[Corank>0,(*If at least one or more eigenvalues of the Hessian are zero, indicated by a non-zero corank, then we have a higher order critical point.*)

If[IntegerQ[MaxDeg]&&MaxDeg>0,(*Sanity check to ensure that we have a legitimate polunomial.*)
PolyDeg=MaxDeg;
,
PolyDeg=ResourceFunction["PolynomialDegree"][Poly,KVec]+1 
(*We add one to the actual degree of the polynomial since the algorithm can narrow down the determinacy to either the actual det or det+1, so that the check for determinacy has to be performed up to at least up one degree higher than the suspected determinacy of the polynomial.*)
];

(*Given polynomial p and variables Subscript[x, 1],Subscript[x, 2],...,Subscript[x, n], we need to compute the first derivatives of the polynomial with respect to the variables:*)
PolyDerivs=First@Last@Reap[
Do[
Sow[D[Poly,Ki]]
,{Ki,KVec}]
];

(*As we attempt to check for determinacy at increasing orders, at each order k, we first accumulate all the monomials of degree k into the list AllMonoms. We then accumulate the basis polynomials obtained by multiplying these monomials by the derivative polynomials in PolyDerivs, into BasisPolys.*)
AllMonoms={};
BasisPolys={};
ZeroList={};

On[LinearSolve::nosol];(*Turning on the notification in case it has been turned off by the user.*)

Do[(*A loop to check for determinacy order by order up to one degree higher than the actual degree of the polynomial (or up to a level specified by the user).*)

(*Our strategy for constructing all the monomials of degree k in n variables is to first compute the distinct partitions of k into n positve integers. Mathematica returns this as a collection of sets of length n. We then repeat the procedure for partitions of k into n-1 positive integers, then n-2 positive integers and so on until partition of k into one positive integer. For those partitions of k into less than n integers (where Mathematica returns sets of length less than n), we append zeros into the set to make up n integers in total. After this, we generate all the permutations of this set, which then gives us all possible ways of writing k = Subscript[i, 1] + Subscript[i, 2] + ... + Subscript[i, n] for integer Subscript[i, j] >= 0. This then allows us to compute all the monomials of degree k in n variables.*)
PostvDistPartns=IntegerPartitions[CurrDeg,{SpaceDim}]; (*Partitions of k into n positive integers.*)

PartnsWithZeros={}; (*Partitions of k into less than n positive integers. Rest of elements are taken to be zeros.*)

Do[
PartnsWithZeros=Join[PartnsWithZeros,(Join[ConstantArray[0,i],#])&/@IntegerPartitions[CurrDeg,{SpaceDim-i}]]
,{i,1,SpaceDim-1}];

Partns=Flatten[Permutations/@Join[PostvDistPartns,PartnsWithZeros],1];

Monoms=First@Last@Reap[
Do[
Sow[Product[KVec[[i]]^CurrPartn[[i]],{i,1,SpaceDim}]]
,{CurrPartn,Partns}]
];

AllMonoms=Join[AllMonoms,Monoms]; (*Appending the newly computed monomials into the list AllMonoms.*)

(*We will be solving the linear system C.B = S, where C is the coefficient matrix of the base polynomials in the basis of monomials(the base polynomials are ontaining by multiplying the first derivatives of the original polynomial by all the monomials. S is a set containing 1's at the position of the monomials of degree k, and zero elsewhere (since the aim is to check whether a linear combination of the base polynomials can yield all the monomials of degree k))*)
OneList=ConstantArray[1,Length[Monoms]];
BMat=Join[ZeroList,OneList];
ZeroList=Join[ZeroList,ConstantArray[0,Length[Monoms]]];

BasisPolys=Join[BasisPolys,
First@Last@Reap[
Do[
Do[
Sow[polyderiv monom//Expand]
,{polyderiv,PolyDerivs}]
,{monom,Monoms}]
]
];

CoeffMatrix=First@Last@Reap[
Do[
CoeffRow=First@Last@Reap[
Do[
Sow[ReplaceAll[Coefficient[poly,monom],Table[Ki->0,{Ki,KVec}]]];
,{monom,AllMonoms}]
];

Sow[CoeffRow]
,{poly,BasisPolys}]
];

SolMat=Quiet[Check[LinearSolve[Transpose[CoeffMatrix],BMat],1]]; (*Solving the linear system. If it has a solution then the polynomial is k-determinate.*)


If[NumberQ[SolMat]&&SolMat==1,
Solved=False;
,
Solved=True;
Determinacy=CurrDeg;

CoeffMatrix=Join[CoeffMatrix,
First@Last@Reap[
Do[
CoeffRow=First@Last@Reap[
Do[
Sow[ReplaceAll[Coefficient[polyderiv,monom],Table[Ki->0,{Ki,KVec}]]];
,{monom,AllMonoms}]
];

Sow[CoeffRow]

,{polyderiv,PolyDerivs}]
]
];

(*Codimension=(Determinacy^2+3Determinacy)/2-MatrixRank[CoeffMatrix]; for the two dimensional case.*)

Codimension=Length[AllMonoms]-MatrixRank[CoeffMatrix];


Print["Polynomial has corank ",Corank,". It is ",CurrDeg,"-determinate. Its determinacy is either ",CurrDeg," or ",CurrDeg-1," and its codimension is ",Codimension,". Comparing with the known list of singularities:"];
Break[]
];

,{CurrDeg,1,PolyDeg}];
,

(*If all elements of the Jacobian are zero and Corank = 0, we have a good old non-singular quadratic form that could be a maximum, minimum or a saddle.*)
Determinacy=2;
Codimension=1;

Solved=True;

];

Switch[{Corank,Codimension,Determinacy},
{0,1,2},Print["Polynomial is equivalent to an ordinary quadratic form with a non-singular Hessian."];,
{1,1,3},Print["Polynomial has a fold (\!\(\*SubscriptBox[\(A\), \(2\)]\)) catastrophe."];,
{1,1,4},Print["Polynomial has a fold (\!\(\*SubscriptBox[\(A\), \(2\)]\)) catastrophe."];,
{1,2,4},Print["Polynomial has a cusp (\!\(\*SubscriptBox[\(A\), \(3\)]\)) catastrophe."];,
{1,2,5},Print["Polynomial has a cusp (\!\(\*SubscriptBox[\(A\), \(3\)]\)) catastrophe."];,
{1,3,5},Print["Polynomial has a swallowtail (\!\(\*SubscriptBox[\(A\), \(4\)]\)) catastrophe."];,
{1,3,6},Print["Polynomial has a swallowtail (\!\(\*SubscriptBox[\(A\), \(4\)]\)) catastrophe."];,
{1,4,6},Print["Polynomial has a butterfly (\!\(\*SubscriptBox[\(A\), \(5\)]\)) catastrophe."];,
{1,4,7},Print["Polynomial has a butterfly (\!\(\*SubscriptBox[\(A\), \(5\)]\)) catastrophe."];,
{1,5,7},Print["Polynomial has a wigwam (\!\(\*SubscriptBox[\(A\), \(6\)]\)) catastrophe."];,
{1,5,8},Print["Polynomial has a wigwam (\!\(\*SubscriptBox[\(A\), \(6\)]\)) catastrophe."];,
{1,6,8},Print["Polynomial has a star (\!\(\*SubscriptBox[\(A\), \(7\)]\)) catastrophe."];,
{1,6,9},Print["Polynomial has a cusp (\!\(\*SubscriptBox[\(A\), \(7\)]\)) catastrophe."];,
{2,3,3},Print["Polynomial has either an elliptic umbilic (\!\(\*SubsuperscriptBox[\(D\), \(4\), \(-\)]\)) (monkey saddle) or a hyperbolic umbilic (\!\(\*SubsuperscriptBox[\(D\), \(4\), \(+\)]\)) catastrophe."];,
{2,3,4},Print["Polynomial has either an elliptic umbilic (\!\(\*SubsuperscriptBox[\(D\), \(4\), \(-\)]\)) (monkey saddle) or a hyperbolic umbilic (\!\(\*SubsuperscriptBox[\(D\), \(4\), \(+\)]\)) catastrophe."];,
{2,4,4},Print["Polynomial has a parabolic umbilic (\!\(\*SubscriptBox[\(D\), \(5\)]\)) catastrophe."];,
{2,4,5},Print["Polynomial has a parabolic umbilic (\!\(\*SubscriptBox[\(D\), \(5\)]\)) catastrophe."];,
{2,5,4},Print["Polynomial has a symbolic umbilic (\!\(\*SubscriptBox[\(E\), \(6\)]\)) catastrophe."];,
{2,5,5},Print["Polynomial has a symbolic umbilic (\!\(\*SubscriptBox[\(E\), \(6\)]\)), or a second elliptic umbilic (\!\(\*SubsuperscriptBox[\(D\), \(6\), \(-\)]\)) or a second hyperbolic umbilic (\!\(\*SubsuperscriptBox[\(D\), \(6\), \(+\)]\)) catastrophe."];,
{2,5,6},Print["Polynomial has either a second elliptic umbilic (\!\(\*SubsuperscriptBox[\(D\), \(6\), \(-\)]\)) or a second hyperbolic umbilic (\!\(\*SubsuperscriptBox[\(D\), \(6\), \(+\)]\)) catastrophe."];,
{2,6,4},Print["Polynomial has a \!\(\*SubscriptBox[\(E\), \(7\)]\) catastrophe."];,
{2,6,5},Print["Polynomial has a \!\(\*SubscriptBox[\(E\), \(7\)]\) catastrophe."];,
{2,6,6},Print["Polynomial has a \!\(\*SubscriptBox[\(D\), \(7\)]\) catastrophe."];,
{2,6,7},Print["Polynomial has a \!\(\*SubscriptBox[\(D\), \(7\)]\) catastrophe."];,
{2,7,5},Print["Polynomial has a \!\(\*SubscriptBox[\(E\), \(8\)]\) catastrophe."];,
{2,7,6},Print["Polynomial has a \!\(\*SubscriptBox[\(E\), \(8\)]\) catastrophe."];,
{2,7,7},Print["Polynomial has a either \!\(\*SubsuperscriptBox[\(D\), \(8\), \(+\)]\) or a \!\(\*SubsuperscriptBox[\(D\), \(8\), \(-\)]\) catastrophe."];,
{2,7,8},Print["Polynomial has a either \!\(\*SubsuperscriptBox[\(D\), \(8\), \(+\)]\) or a \!\(\*SubsuperscriptBox[\(D\), \(8\), \(-\)]\) catastrophe."];,
_,
If[Solved==True,
If[Corank==2&&Codimension==8&&(Determinacy==4||Determinacy==5)&&SpaceDim==2,
Print["Polynomial has a catastrophe belonging to the \!\(\*SubscriptBox[\(X\), \(9\)]\) singularity class."];
,
Print["Point singularity has a finite determinacy but cannot identified by name!"]
],
Print["The polynomial appears to host a singularity that does not have a finite determinacy."];
Determinacy=\[Infinity];
Codimension=\[Infinity];
]
]
,

Solved=True;
Determinacy=1;
Codimension=0;
Print["The given polynomial does not have a critical point at the origin. It is smoothly equivalent to a linear form, and has determinacy 1."];
];

{Corank,Codimension,Determinacy}

]


(* ::Section:: *)
(*Computing the series expansion of the eigenvalues of a matrix*)


(*This function is still in the works. It needs to be expanded to carefully incorporate the extension to non-trivially degenerate bands. In this version, it can handle situations where the degenerate bands have same Taylor expansion up to the degree we are interested in.*)

Options[TaylorExpandBand]={"Dimension"->-1,"TuningDegree"->0,"TuningSymbol"->"\[Delta]t","Resolution"->10^-8,"SampleDirection"->-1,"Messages"->True}

TaylorExpandBand[Hamltn_,k0_,BandNo_,TaylorDegreeIn_,OptionsPattern[]]:=
Module[{SpaceDimIn,HoppingParamDegree,TuningParamSymbol,Resln,DegenDirn,YesMessage,Eigs,EigVals,EigVecs,HamltnDim,HamltnDerivs,ValidHamltn,Broken,DegenBreakingOrder,DiagonaliseOrder,HMatrix,EMatrix,ProjEMatrix,ProjEigVecs,DegenNK,DerivMatrix,DerivMatrixConj,EDerivs,ThMatrix,TempSet,TempSetCopy,CurrBnd,NewSubSet,DiagonalComponent,ComplementSet,TaylorPoly,TaylorPolyArray,HoppingDim,TaylorDegree,kVec,tVec,ktVec,gn,fn,\[Lambda]Vec,Degen,NoDegen,NonDegen,NoNonDegen,SpaceDim,HoppingParamDegreeFinal,SampleDirn,Diagonalised,ExtraTerm,ExtraTermRefined,EDerivsPolished,RefSet,EigTrfmn},

SpaceDimIn=OptionValue["Dimension"];
HoppingParamDegree=OptionValue["TuningDegree"];
TuningParamSymbol=OptionValue["TuningSymbol"];
Resln=OptionValue["Resolution"];
DegenDirn=OptionValue["SampleDirection"];
YesMessage=OptionValue["Messages"];

If[Length[k0]>0&&AllTrue[k0,NumberQ[#]&]&&IntegerQ[BandNo]&&BandNo>0&&IntegerQ[TaylorDegreeIn]&&TaylorDegreeIn>0&&IntegerQ[HoppingParamDegree]&&HoppingParamDegree>=0&&IntegerQ[SpaceDimIn]&&NumberQ[Resln]
(*Some sanity checks to ensure that valid inputs have been fed in by the user. The point k0 is a vector containing both the k-point and the values of the tuning parameters, about which we will Taylor expand. i.e k0 = (Subscript[k, 1],Subscript[k, 2],...,Subscript[k, n],Subscript[t, 1],Subscript[t, 2],...,Subscript[t, m]).*)
,
If[SpaceDimIn==-1||SpaceDimIn>Length[k0],
(* If either the k-space dimension is not provided (in which case it defaults to -1), or if the provided value does not make sense when compared to dimension of k0, then we interpret the dimension of the k-point/vector of interest k0 to be the dimension of the full k-space (i.e there are no parameters t provided). *)

SpaceDim=Length[k0]; 

HoppingDim=0; (* As a reflection of the statement above, none of the components in k0 represent hopping parameters, so that we set the dimension of the hopping sub-vector to zero. *)
HoppingParamDegreeFinal=0;

,
SpaceDim=SpaceDimIn;
HoppingDim=Length[k0]-SpaceDimIn ;
(* The number of hopping parameters is equal to dim k0 - dim of k-space. Recall that the input k0 contains both a k-part and a t-part.*)
HoppingParamDegreeFinal=HoppingParamDegree;
];

(*We now define symbolic k and t vectors which will be used for the series expansion.*)
kVec=Table[Symbol["p"<>ToString[i]],{i,1,Length[k0]}]; 

If[StringQ[TuningParamSymbol], (*If one valid string symbol is provided, we use it with the appropriate subscript for all the parameters.*)
If[HoppingDim==1,
tVec={Symbol[StringDelete[TuningParamSymbol," "]]}; (*obviously we have to remove any whitespaces!*)
,
tVec=Table[Symbol[StringDelete[TuningParamSymbol," "]<>ToString[i]],{i,1,HoppingDim}];
];
,
If[Length[TuningParamSymbol]==HoppingDim&&AllTrue[TuningParamSymbol,StringQ[#]&],
tVec=Symbol/@((StringDelete[#," "])&/@TuningParamSymbol); (*Deleting whitespace and then coverting them to symbols.*)
,
tVec=Table[Symbol["\[Delta]t"<>ToString[i]],{i,1,HoppingDim}];
]
];

ktVec=Join[Table[Symbol["p"<>ToString[i]],{i,1,SpaceDim}],tVec];
(*Print["Full (k,t) symbol set = ",ktVec,"."];*)

\[Lambda]Vec={Symbol["\[Lambda]k"],Symbol["\[Gamma]t"]};

If[YesMessage==True,
Print["\!\(\*
StyleBox[\"k\",\nFontWeight->\"Bold\",\nFontSlant->\"Italic\"]\)-space dimension = ",SpaceDim,"."];
Print["Parameter space (\!\(\*
StyleBox[\"t\",\nFontWeight->\"Bold\",\nFontSlant->\"Italic\"]\)-space) dimension = ",HoppingDim,"."];
Print["We shall attempt to series expand the given Hamiltonian at ",ktVec," = ",k0," to degree ",TaylorDegreeIn," and ",HoppingParamDegree," in \!\(\*
StyleBox[\"k\",\nFontWeight->\"Bold\",\nFontSlant->\"Italic\"]\) and \!\(\*
StyleBox[\"t\",\nFontWeight->\"Bold\",\nFontSlant->\"Italic\"]\) respectively."];
];

(*In the case of degeneracies, we will attempt to diagonalise a matrix that is symbolic. To circumvent issues, we will attempt to diagonalise the symbolic matrix at some sample numeric kVec (i.e sample direction) and check if the eigensystem manages to diagonalise the symbolic matrix as well.*)
If[Length[DegenDirn]==Length[k0]&&AllTrue[DegenDirn,NumberQ[#]&],
If[YesMessage==True,
Print["A valid sample (\!\(\*
StyleBox[\"k\",\nFontWeight->\"Bold\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontWeight->\"Plain\"]\)\!\(\*
StyleBox[\"t\",\nFontWeight->\"Bold\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\")\",\nFontWeight->\"Plain\"]\)-vector ",DegenDirn," seems to have been provided. If we are lucky, we will be able to diagonalise the symbolic matrix ",Subsuperscript["\[ScriptCapitalE]","N","(\[Kappa])"]," by first evaluating it at this sample vector and then numerically diagonalising it."];
];
SampleDirn=DegenDirn;
,
SampleDirn=ConstantArray[0,Length[k0]];
Do[
SampleDirn+=RandomReal[] UnitVector[Length[k0],i];
(*We are essentially populating the vector components by pseudo random numbers between 0 and 1.*)
,{i,1,Length[k0]}];

If[YesMessage==True,
Print["Invalid sample (\!\(\*
StyleBox[\"k\",\nFontWeight->\"Bold\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontWeight->\"Plain\"]\)\!\(\*
StyleBox[\"t\",\nFontWeight->\"Bold\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\")\",\nFontWeight->\"Plain\"]\)-vector provided for numerical diagonailsation of symbolic matrix. Using the randomly generated ",SampleDirn," instead."];
];
];

TaylorDegree=TaylorDegreeIn+HoppingParamDegreeFinal;
(*This is needed since at the k order equal to TaylorDegree, we would not have any terms involving the the tuning parameters t's if we did not expand to order TaylorDegree + HoppingParamDegree to begin with.*)

		
(*Checking if a non-symbolic Hamiltonian has been fed:*)
ValidHamltn=True;
Do[
If[!AllTrue[HRow,NumberQ[#]&],
If[YesMessage==True,
Print[Style["The given Hamiltonian, evaluated at (\!\(\*
StyleBox[\"k\",\nFontWeight->\"Bold\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"t\",\nFontWeight->\"Bold\",\nFontSlant->\"Italic\"]\)) = ",k0," has symbolic terms. We are unable to compute the series expansion!",Red]];
];

ValidHamltn=False;
]
,{HRow,N[Hamltn[k0]]}];

If[!ValidHamltn,
Return[]
];

(*Checking if the Hamiltonian matrix at k0 is Hermitian. Skip if not needed*)
(*
If[!HermitianMatrixQ[Chop[N[Hamltn[k0]],Resln]],
Print[Style["The given Hamiltonian is not hermitian. Series expansion can not be evaluated!",Red]];
Return[]
];
*)

(* Computing the eigenvalues and eigenvectors at k0. We use 1/2(H+H^\[Dagger]) to ensure Hermiticity.*)
Eigs=Eigensystem[0.5 (N[Hamltn[k0]]+Transpose[Conjugate[N[Hamltn[k0]]]]),Method->"Direct"]; 

EigVals=Eigs[[1]];
EigVecs=Eigs[[2]];

(* For the calculations that will follow, we need the dimension of the Hamiltonian matrix, which is the same as the number of bands: *)
HamltnDim=Length[EigVals];

Eigs=SortBy[Table[Append[{EigVals[[i]]},EigVecs[[i]]],{i,1,HamltnDim}],First]; (* Sorting the eigensystem by the eigenvalues, rather than by their magnitudes, which is what Mathematica does by default. *)

EigVals=Eigs[[All,1]];
EigVecs=Eigs[[All,2]];
(* Separating the eigenvalues and eigenvectors in to separate lists *)

If[YesMessage==True,
Print["Precision of eigensystem solution at ",k0,", evaluated as Norm[\!\(\*SuperscriptBox[\(EigVecs\), \(*\)]\).H.\!\(\*SuperscriptBox[\(EigVecs\), \(T\)]\)-DiagonalMatrix[EigVals]] = ",Norm[Conjugate[EigVecs] . (0.5 (N[Hamltn[k0]]+Transpose[Conjugate[N[Hamltn[k0]]]])) . Transpose[EigVecs]- DiagonalMatrix[EigVals]],"."];

Print["The band eigenvalues: ",EigVals,"."];
];

(* We now attempt to find the subset of eigenvectors that are degenerate to the band of interest, which is denoted by BandNo. Obviously, the degeneracy is determined by the resolution adopted by us. *)
Degen=First@Last@Reap[Do[
If[Abs[EigVals[[i]] - EigVals[[BandNo]]]<=Resln,Sow[i]];
,{i,1,HamltnDim}]];

NoDegen=Length[Degen];

(* Having found the subset of bands dengenerate to the band of interest, we now find the subset of bands that are not degenerate to it: *)
NonDegen=Complement[Table[i,{i,1,HamltnDim}],Degen];

DegenNK={{Table[j,{j,1,NoDegen}]}}; 
(*We use a relative indexing starting from 1 for the D_N^K sets.
The first element is D_N^0 which is the set of all the degenerate bands, 
which we will index starting from 1. For example {{1,2,3,4,5,6}}. 
The reason for double {{}} will become clear below.

Subsequent elements will contain D_N^K with the relevant groupings, for example
{{1,3,4},{2,6},{5}} etc
*)

ThMatrix={}; (*The \[CapitalTheta]_N^{-K} matrices.*)

(* Number of bands that are not degenerate to BandNo: *)
NoNonDegen=Length[NonDegen];

(* Some sanity checks...can be removed if needed. *)
If[YesMessage==True,
Print["Degenerate band(s) to be series expanded: ",Degen,"."];
Print["Rest of the band(s): ",NonDegen,"."];
];

DegenBreakingOrder=-1; (*This is the order k at which \[PartialD]^kH projected to the degenerate subspace is not diagonal.*)

(* Lists of matrices that we will be using for the Taylor expansion calculation: *)
HamltnDerivs=ConstantArray[0,TaylorDegree];
HMatrix=ConstantArray[0,{TaylorDegree,Length[Degen]}];
EMatrix=ConstantArray[0,TaylorDegree];
EDerivs=ConstantArray[0,{TaylorDegree+1,Length[Degen]}];
EDerivs[[1]]=Part[EigVals,Degen];
(* EDerivs contains the polynomials \[PartialD]^lSubscript[\[Epsilon], n] that make up the Taylor expansion, stored order by order. It is a multidimensional list since we will be computing the Taylor expansion of each of the degenerate bands. If the bands have the same Taylor expansion to the given order at k0, then this is a redundancy. This means that their Taylor expansions will be identical. *)

(* The following is a function that will be needed in the calculation to follow. Basically it computes the quantity \[LeftAngleBracket]Subscript[m, i]|\[PartialD]^(Subscript[k, i]-Subscript[k, i+1])H|Subscript[m, i+1]\[RightAngleBracket]-\[PartialD]^(Subscript[k, i]-Subscript[k, i+1])Subscript[\[Epsilon], n].Subscript[\[Delta], Subscript[m, i],Subscript[m, i+1]] :*)
gn[mi_,ki_,mip1_,kip1_,ni_]:=If[1<=mi<=HamltnDim&&1<=mip1<=HamltnDim&&1<=kip1<ki&&MemberQ[Degen,ni],HamltnDerivs[[ki-kip1(*+1*),mi,mip1]]-KroneckerDelta[mi,mip1]EDerivs[[ki-kip1+1,First@First@Position[Degen,ni]]],0];
(* The +1 is needed in the array indices since Mathematica indexing begins from 1 rather than 0. The first index then represents 0 in the conventional sense. Also we don't have to worry about the case when kip1 = ki, since it will not occur. The function appears inside sums where kip1 runs from 1 to ki-1. This is the reason why we check kip1 < ki and not kip1 <= ki. Lastly, we use First@First@Position[Degen,n] instead of n directly since the EDerivs array is indexed by 1,...,Length[Degen] and n instead represents the index of the band. *)

fn[ki_,mi_,li_,ni_]:=If[MemberQ[Degen,ni]&&MemberQ[NonDegen,mi]&&1<=li<=HamltnDim&&ki>=1,1/(EigVals[[ni]]-EigVals[[mi]]) (HamltnDerivs[[ki(*+1*),mi,li]]+Sum[Binomial[ki,kip1]Sum[gn[mi,ki,mip1,kip1,ni]fn[kip1,mip1,li,ni],{mip1,NonDegen}],{kip1,1,ki-1}]),0];


(* We now compute, order by order, the set of monomials that make up the Taylor expansion of the Hamiltonian, all the way up to the order specified by TaylorDegree. For the convenience of the forthcoming calculations, it makes sense to transform all these into the basis of the eigenvectors at k0: *)

(*Print["Hamiltonian matrix in the eigenbasis ",Chop[Conjugate[EigVecs] . (1/2(N[Hamltn[k0]]+Transpose[Conjugate[N[Hamltn[k0]]]])). Transpose[EigVecs],Resln]];*)

Broken=False;
Diagonalised=False;
ExtraTerm=ConstantArray[0,Length[Degen]];


Do[
If[YesMessage==True,
Print[Superscript["D","("<>ToString[(i-1)]<>")"]," = ",(Degen[[#]])&/@DegenNK[[i]]];
];

DerivMatrix=Simplify[N[D[Hamltn[l*kVec+k0//N],{l,i}]/.{l->0}]];

DerivMatrixConj=ComplexExpand[Conjugate[Transpose[DerivMatrix]]];

(*The final form of DerivMatrix at any order k should be quite simple since each of its terms contains only sum of monomial terms of degree k. If we don't apply simplify in the lines above, Mathematica will instead have a messy soup that will complicate the calculations below.*)

HamltnDerivs[[i]]=Simplify[Conjugate[EigVecs] . (0.5(DerivMatrix+DerivMatrixConj) ) . Transpose[EigVecs]];
(* Normally we would do basis transformations of a matrix M by Evec^\[Dagger].M.Evec. However, since Mathematica stores the eigenvectors in the form of rows, we have to instead do (Evec)^*.M.Evec^T. *)

(*We shall be checking if ProjHamltnDerivs, which is the projection of the derivative(s) of the Hamiltonian to the degenerate subspace contains only diagonal matrices.*)
Do[
HMatrix[[i,DegenCurr]]=First@Last@Reap[
Do[
Sow[
First@Last@Reap[
Do[
Sow[Simplify[HamltnDerivs[[i,ind1,ind2]]+Sum[Binomial[i,k1]Sum[HamltnDerivs[[i-k1,ind1,m1]]fn[k1,m1,ind2,Degen[[DegenCurr]]],{m1,NonDegen}],{k1,1,i-1}]]];
,{ind2,Degen}]
]
]
,{ind1,Degen}]
];

HMatrix[[i,DegenCurr]]=0.5(HMatrix[[i,DegenCurr]]+ComplexExpand[ConjugateTranspose[HMatrix[[i,DegenCurr]]]])//Simplify;

,{DegenCurr,1,Length[Degen]}];

(*
If[i<=2,
(*Print["H deriv = ",Chop[HamltnDerivs[[i]],Resln]];*)
Print["HMatrix = ",Chop[HMatrix[[i]],Resln]];
];
*)

(*The \[ScriptCapitalE] matrix is defined now: *)
If[i!=4,
EMatrix[[i]]={};
Do[
AppendTo[EMatrix[[i]],HMatrix[[i,DegenNK[[i,j,1]]]]]; 
(*Print["Appended!"];*)
(*Like HMatrix, EMatrix is a list of a list of matrices at each order.*)
,{j,1,Length[DegenNK[[i]]]}];

,

EMatrix[[i]]={};
Do[
If[YesMessage==True,
Print["Computing the corrections to ",Subsuperscript["\[ScriptCapitalE]",Degen[[DegenNK[[i,j,1]]]],"("<>ToString[i]<>")"],"..."];
];

ExtraTerm=4HMatrix[[3,DegenNK[[i,j,1]]]] . ThMatrix[[2,Position[DegenNK[[i-2+1]],DegenNK[[i,j,1]]][[1,1]]]] . HMatrix[[3,DegenNK[[i,j,1]]]];
(*To identify the correct \[CapitalTheta]^{(2)} matrix to use, we have to pick a sample element from the j^th D_N^{i-1} subset = DegenNK[[i,j,1]].
We have to identify its D_N^2 subset, that is we have to compute its position in DegenNK[[i-2+1]].
*)

ExtraTermRefined=First@Last@Reap[
Do[
Sow[
First@Last@Reap[
Do[
Sow[Chop[Expand[PolynomialReduce[Numerator[Together[ExtraTerm[[ind1,ind2]]]],Denominator[Together[ExtraTerm[[ind1,ind2]]]],kVec][[1,1]]],Resln]];
(*
If[!NumberQ[ExtraTerm[[ind1,ind2]]],
(*Print[Numerator[Chop[ExtraTerm[[ind1,ind2]],Resln]]];
Print[Denominator[Chop[ExtraTerm[[ind1,ind2]],Resln]]];*)
Print[Chop[Expand[PolynomialReduce[Numerator[Together[ExtraTerm[[ind1,ind2]]]],Denominator[Together[ExtraTerm[[ind1,ind2]]]],kVec][[1,1]]],Resln]]
]*)
,{ind2,1,NoDegen}]
]
]
,{ind1,1,NoDegen}]
];

(*Print[Chop[ExtraTermRefined,Resln]//MatrixForm];*)

AppendTo[EMatrix[[i]],(HMatrix[[i,DegenNK[[i,j,1]]]]+ExtraTermRefined)];

If[YesMessage==True,
Print["Done!"]
];
,{j,1,Length[DegenNK[[i]]]}];
];

(*Now that the \[ScriptCapitalE] matrices are defined, we have to diagonalise them in the appropriate subspaces, i.e the D_N^{i-1} = DegenNK[[i]].*)
(*We first define a matrix that will later be used to basis-transform the HMatrices, EMatrices and HamltnDeriv matrices.*)
EigTrfmn=ConstantArray[0,{HamltnDim,HamltnDim}];

If[NoNonDegen>0,
Part[EigTrfmn,NonDegen,NonDegen]=IdentityMatrix[NoNonDegen];
];
(*As for the matrix elements between the degenerate bands, they shall be filled based on the diagonalisations performed below.*)

Do[
(*We are iterating through different D_N^{(i-1)}, that is DegenNK[[i,j]]'s.
We project \[ScriptCapitalE]_N^{i} to this subspace and diagonalise it at the sample k vector.
*)
ProjEMatrix=Part[EMatrix[[i,j]],DegenNK[[i,j]],DegenNK[[i,j]]];
ProjEigVecs=Eigenvectors[N[ReplaceAll[ProjEMatrix,Table[kVec[[l]]->SampleDirn[[l]],{l,Length[kVec]}]]//Simplify]];

ProjEMatrix=Simplify[Conjugate[ProjEigVecs] . ProjEMatrix . Transpose[ProjEigVecs]];

(*We will now check if this is diagonal*)
Diagonalised=True;(*Provisionally*)

Do[
Do[
If[!NumberQ[Chop[ProjEMatrix[[m,n]],Resln]],
	If[YesMessage==True,
	Print["Non-zero off diagonal term equal to ",Chop[ProjEMatrix[[m,n]],Resln]," encountered at order ",i,". We were unable to diagonalise the symbolic matrix by diagonalising it along the ",ktVec," = ",SampleDirn," direction. Try feeding in a different direction/vector and recalling the function afresh."];
	Print["The degenerate subset is ",DegenNK[[i,j]]," and the relative location of matrix element is ",m,", ",n,"."];
	];
	Diagonalised=False;
]
,{n,m+1,Length[ProjEMatrix]}]
,{m,1,Length[ProjEMatrix]}];

If[Diagonalised,

If[YesMessage==True,
Print["The symbolic, projected \"derivative\" matrix ",Subsuperscript["\[ScriptCapitalE]",Degen[[DegenNK[[i,j,1]]]],"("<>ToString[i]<>")"]," was successfully diagonalised."];
];
(*Print["The eigenvalues of the symbolic matrix are ",Chop[Diagonal[ProjHamltnDerivsCurr],Resln]];*)
(*Now that we have successfully diagonalised the symbolic matrix, we have to 'rotate' the eigenvectors in the degenerate subspace.*)

(*Since E^\[Dagger].(\[CapitalDelta]^lH).E is diagonal, the new eigenvectors are given by (Subscript[v, 1], Subscript[v, 2],..., Subscript[v, n]).E, where Subscript[v, i] are the original degenerate eigenvectors. These are in the column vector system. In case of row based system as in Mathematica, we would instead have E.(Subscript[v, 1], Subscript[v, 2],..., Subscript[v, n])^T.*)
RefSet=(Degen[[#]])&/@DegenNK[[i,j]];
(*Print[DegenNK[[i,j]],"->",RefSet];*)
Part[EigVecs,RefSet]=ProjEigVecs . Part[EigVecs,RefSet];
(*We are using (Degen[[#]])&/@DegenNK[[i,j]] to find the correct band indices since we have been using a relative index system within the degenerate subset.*)

Part[EigTrfmn,RefSet,RefSet]=Transpose[ProjEigVecs];

If[YesMessage==True,
Print["Precision of the \"rotated\" eigensystem solution = ",Norm[Conjugate[EigVecs] . (0.5 (N[Hamltn[k0]]+Transpose[Conjugate[N[Hamltn[k0]]]])) . Transpose[EigVecs]- DiagonalMatrix[EigVals]],"."];
(*Print["The diagonalised matrix: ",Chop[ProjEMatrix,Resln]];*)
];

Part[EDerivs,i+1,DegenNK[[i,j]]]=Diagonal[ProjEMatrix]

(*Part[EMatrix[[i,j]],DegenNK[[i,j]],DegenNK[[i,j]]]=ProjEMatrix*)
(*We still have to recompute the derivative matrix and HMatrix at order i since the basis has changed. But we will do it once all these basis change computations are done and we have the final eigenbasis at order i.*)
];

,{j,1,Length[DegenNK[[i]]]}];

(*
If[i<=2,
Print[Chop[MatrixForm[EigTrfmn],Resln]]
];*)

(*Since we diagonalised the \[ScriptCapitalE] matrices in the degenerate subspaces, we have to recompute the derivative matrices and the HMatrices again all the way from beginning to now.*)
Do[
HamltnDerivs[[i2]]=Simplify[ConjugateTranspose[EigTrfmn] . HamltnDerivs[[i2]] . EigTrfmn];

Do[
HMatrix[[i2,j2]]=Simplify[ConjugateTranspose[Part[EigTrfmn,Degen,Degen]] . HMatrix[[i2,j2]] . Part[EigTrfmn,Degen,Degen]];
,{j2,1,Length[Degen]}];

Do[
EMatrix[[i2,j2]]=Simplify[ConjugateTranspose[Part[EigTrfmn,Degen,Degen]] . EMatrix[[i2,j2]] . Part[EigTrfmn,Degen,Degen]];
,{j2,1,Length[DegenNK[[i2]]]}]

,{i2,1,i}];

(*HMatrix[[i]]=0.5(HMatrix[[i]]+ComplexExpand[ConjugateTranspose[HMatrix[[i]]]])//Simplify;*)
(*
If[i<=2,
Print["HMatrix after rotation: ",Chop[HMatrix[[i]],Resln]]
];*)

(*Finally we will construct D_N^i which we denote by DegenNK[[i+1]]*)
AppendTo[DegenNK,{}]; (*This will be DegenNK[[i+1]], which we shall populate below.*)
AppendTo[ThMatrix,{}];

Do[
TempSet=DegenNK[[i,j]]; (*Recall that like EDerivs, DegenNK also has a zeroth order element which in Mathematica notation is DegenNK[[1]].*)

While[Length[TempSet]>0,(*If TempSet has become empty, we have iterated through all its elements and allocated them to the appropriate i'th-degenerate subsets.*)
CurrBnd=TempSet[[1]];
TempSet=Delete[TempSet,1];
TempSetCopy=TempSet;
NewSubSet={CurrBnd};(*NewSubSet contains all the elements in TempSet which are i'th-degenerate to CurrBnd*)

Do[
If[Chop[(EDerivs[[i+1,CurrBnd]]-EDerivs[[i+1,TempSetElmt]])//Simplify,Resln]==0, (*i+1 needed in the index since EDerivs[[1]] cotains the zeroth derivatives, that is the eigenvalues themselves.*)
TempSetCopy=Delete[TempSetCopy,Position[TempSetCopy,TempSetElmt][[1,1]]];

AppendTo[NewSubSet,TempSetElmt];
]
,{TempSetElmt,TempSet}];

TempSet=TempSetCopy;
AppendTo[DegenNK[[i+1]],NewSubSet];

ComplementSet=Complement[DegenNK[[i,j]],NewSubSet];

DiagonalComponent=ConstantArray[0,NoDegen];
Do[
DiagonalComponent[[Elmt]]=1/(Binomial[i+1,i](EDerivs[[i+1,NewSubSet[[1]]]]-EDerivs[[i+1,Elmt]]));
,{Elmt,ComplementSet}];

AppendTo[ThMatrix[[i]],DiagonalMatrix[DiagonalComponent]];
];

,{j,1,Length[DegenNK[[i]]]}]; (*iterating through the elements of DegenNK[[i-1]].*)

(*Print[DegenNK[[i+1]]];*)
,{i,1,TaylorDegree}];



If[Diagonalised,

(*Having defined the functions necessary for our calculations, we go ahead and compute the derivative polynomials of the dispersion, denoted by \[PartialD]^lSubscript[\[Epsilon], n], order by order:*)

Do[
Do[
(*EDerivs[[M+1,First@First@Position[Degen,n]]]==Simplify[HamltnDerivs[[M+1,n,n]]+Sum[Binomial[M,k1]Sum[HamltnDerivs[[M-k1+1,n,m1]]fn[k1,m1,n],{m1,NonDegen}],{k1,1,M-1}]];*)

(*We now ensure that no numerically small imaginary terms come in:*)
EDerivs[[M+1,n]]=0.5(EDerivs[[M+1,n]]+ComplexExpand[Conjugate[EDerivs[[M+1,n]]]]);

,{M,1,TaylorDegree}];
,{n,NoDegen}];


(*Before processing the Taylor expansion and assembling it into a format we desire, we first construct it using the polynomials calculated above.*)
TaylorPoly=ConstantArray[0,Length[Degen]];

Do[
Do[
EDerivsPolished=Chop[ReplaceAll[(EDerivs[[M+1,n]]/Factorial[M] (*+ExtraTerm[[n,1,1]]/Factorial[DegenBreakingOrder+2]*)),Join[Table[kVec[[i]]->\[Lambda]Vec[[1]] kVec[[i]],{i,1,SpaceDim}],Table[kVec[[j]]->\[Lambda]Vec[[2]] tVec[[j-SpaceDim]],{j,SpaceDim+1,Length[k0]}]]],Resln];
TaylorPoly[[n]]=TaylorPoly[[n]]+Normal[Series[EDerivsPolished,{\[Lambda]Vec[[2]],0,HoppingParamDegree}]]/.{\[Lambda]Vec[[2]]->1};

(*
If[M==4,
Print["Original term: ",Chop[EDerivs[[M+1,n]],Resln]];
Print["Replacement rules ",Join[Table[kVec[[i]]->\[Lambda]Vec[[1]] kVec[[i]],{i,1,SpaceDim}],Table[kVec[[j]]->\[Lambda]Vec[[2]] tVec[[j-SpaceDim]],{j,SpaceDim+1,Length[k0]}]]];
Print["After replacement: ",EDerivsPolished];
Print[Chop[Normal[Series[EDerivsPolished,{\[Lambda]Vec[[2]],0,HoppingParamDegree}]]/.{\[Lambda]Vec[[2]]->1},Resln]];
]*)

,{M,0,TaylorDegree}];

,{n,1,Length[Degen]}];

(*Now we process the Taylor expansion into a nice and convenient form:*)
TaylorPolyArray=ConstantArray[0,{Length[Degen],TaylorDegreeIn+1}];

Do[
Do[
TaylorPolyArray[[n,i+1]]=1/Factorial[i] Chop[D[TaylorPoly[[n]],{\[Lambda]Vec[[1]],i}]/.{\[Lambda]Vec[[1]]->0},Resln]//FullSimplify;
,{i,0,TaylorDegreeIn}]
,{n,1,Length[Degen]}];

(*Time to return the Taylor expansion, stored as a list, order by order in k:*)
Chop[TaylorPolyArray,Resln]
,
Print[Style["Some issue encountered at order ",Red],Style[ToString[DegenBreakingOrder],Red],Style[". Algorithm did not evaluate successfully!",Red]];
]
,
If[YesMessage==True,
Print[Style["Invalid inputs!",Red]]
]
]
]



(* ::Section::Closed:: *)
(*Binary Search through a sorted list indexed by vectors*)


BinSearch[LocsList_,FindElmt_]:=
Module[
{Found,Done,Lower,Higher,Beg,Final,ElmtPosition,Mid},

Found=False; (*To check if the lattice point is in the list*)
Done=False; (*The search is completed either if the lattice point has been identified or if the apropriate location to insert has been found*)

(*The following are variables used to find the appropriate position for insertion in the event that FindElmt is not already in LocsList*)
Lower=False;
Higher=False;

Beg=1;
Final=Length[LocsList];

(*LocsList is assumed to be already sorted. We check if FindElmt is lower in order than the first element or higher in order than the last element*)
If[Order[LocsList[[1]],FindElmt]==-1,
Done=True;
ElmtPosition=1;
,
If[Order[LocsList[[-1]],FindElmt]==1,
Done=True;
ElmtPosition=Length[LocsList]+1;
]
];


While[
!Done
,
If[Final>=Beg,
Mid=Floor[(Beg+Final)/2];
If[Order[LocsList[[Mid]],FindElmt]==0,
Found=True;
Done=True;
ElmtPosition=Mid;
,
If[Order[LocsList[[Mid]],FindElmt]==1,
Lower=True;
Higher=False;
Beg=Mid+1;
,
Lower=False;
Higher=True;
Final=Mid-1;
]
]
,
(*If Final < Beg, then the element was not found. In particular, in the previous iteration, we had Beg = Final = Mid. Either Beg was increased by 1, or Final was decreased by 1. If Lower is True, then the element at Mid is lower in Order than FindElmt and the appropriate position for insertion is Mid+1. If Higher is true instead, then the element at Mid is higher in order than FindElmt and the appropriate position is Mid.*)
Done=True;

If[Lower,
ElmtPosition=Mid+1;
,
ElmtPosition=Mid;
]
]
];

{ElmtPosition,Found}
]


(* ::Section::Closed:: *)
(*Checking the validity of a symmetrization scheme*)


CheckSymmetrization[SymmList_,BasisVectors_,HamltnDim_,SpaceDim_,YesMessage_,Resln_]:=
Module[{ValidScheme,SymmNewBasis},

ValidScheme=True;

Do[
If[!AllTrue[Flatten[Symm]//N,NumberQ[#]&]||Dimensions[Symm[[1]]]!={SpaceDim,SpaceDim}||Length[Symm[[2]]]!=SpaceDim||Dimensions[Symm[[3]]]!={HamltnDim,HamltnDim},
ValidScheme=False;

If[YesMessage==True,
Print[Style["Each element of the symmetry list must contain an orthogonal matrix acting on real space, a translation vector and a unitary matrix acting on the basis orbitals!",Red]];
];

Break[];
,

SymmNewBasis=Transpose[BasisVectors] . Symm[[1]] . Inverse[Transpose[BasisVectors]];
(*This ensures that the real space symmetry is recast as an orthogonal transformation*)

If[Rationalize[Chop[Transpose[SymmNewBasis] . SymmNewBasis,Resln]]!=IdentityMatrix[SpaceDim]||Rationalize[Chop[ConjugateTranspose[Symm[[3]]] . Symm[[3]],Resln]]!=IdentityMatrix[HamltnDim],
ValidScheme=False;

If[YesMessage==True,
Print[Style["Each element of the symmetry list must contain an orthogonal matrix acting on real space, a translation vector and a unitary matrix acting on the basis orbitals!",Red]];
];

Break[];
]
]
,{Symm,SymmList}];

ValidScheme
]


(* ::Section::Closed:: *)
(*Loading and symmetrizing a tight binding model stored in CSV or Wannier90 formats*)


Options[LoadTBM]={"FileType"->"CSV","LatticeVectors"->-1,"Dimension"->-1,"SpaceDimension"->3,"Messages"->True,"OrbitalCentres"->-1,"Symmetries"->-1,"Resolution"->10^-15,"WannierFileOffset"->4,"SublatticePhase"->False}

(*This funtion loads a tight binding model, stored in a csv or space sparated data file into a Hamiltonian.*)
LoadTBM[DataFile_,OptionsPattern[]]:=Module[
{Strm,SublattPh,CurrLine,TempLine,NoLatticePts,ReadAttempts,LinesRead,DataLength,LineOffset1,MaxLineOffset,BasisVectors,LattPtsList,LattPtLocn,NewLattPt,CurrLattPt,Located,TBMat,TBMatEls,LattVect,indx,NoLines,NewElmt,NeedsSymmetrization,SymmList,AtomLocn,Normalisation,TBMatSymm,TBMatExtras,TBData,Hamltn,kVec,MaxIndex,RowLength,HamltnDim,SpaceDim,YesMessage,Resln,NormalisationExtras},

SpaceDim=OptionValue["SpaceDimension"];
(*YesMessage=OptionValue["Messages"];*)
If[BooleanQ[OptionValue["Messages"]],
YesMessage=OptionValue["Messages"];
,
YesMessage=True;
];

HamltnDim=OptionValue["Dimension"];
AtomLocn=OptionValue["OrbitalCentres"];
SymmList=OptionValue["Symmetries"];
Resln=OptionValue["Resolution"];
MaxLineOffset=OptionValue["WannierFileOffset"];

(*Print[OptionValue["LatticeVectors"]];*)

If[!ListQ[OptionValue["LatticeVectors"]]||!AllTrue[Flatten[OptionValue["LatticeVectors"]]//N,NumberQ[#]&]||!AllTrue[OptionValue["LatticeVectors"],Length[#]==SpaceDim&],
(*If no basis vectors have been provided or if invalid ones are provided, we assume that it is square or cubic lattice*)
If[YesMessage==True,
Print["Invalid or no lattice vectors provided! Using unit vectors."];
];

BasisVectors=IdentityMatrix[SpaceDim];
,
If[YesMessage==True,
Print["Valid lattice vectors provided."];
];

BasisVectors=OptionValue["LatticeVectors"];
];

(*Print[BasisVectors];*)

If[OptionValue["FileType"]=="CSV",

TBData=Import[DataFile,"CSV"];

,
If[OptionValue["FileType"]=="Wannier90",

Strm=OpenRead[DataFile];

CurrLine=0;

(*Skipping the first few lines that just contain the record of when the file was written.*)
Do[
(*Print["Strange stuff ",StringSplit[ReadLine[Strm]][[2]]];*)
TempLine=ToExpression[StringSplit[ReadLine[Strm]][[1]]];
CurrLine++;

If[IntegerQ[TempLine],
HamltnDim=TempLine;
Break[];
]
,{i,1,MaxLineOffset}];

(*Reading the number of lattice points over which the hopping terms are specified.*)
NoLatticePts=ToExpression[ReadLine[Strm]];
CurrLine++;

(*Skipping more lines which have the form 1 1 1 1 1 1 ....*)
If[IntegerQ[NoLatticePts/15],
LineOffset1=IntegerPart[NoLatticePts/15];
,
LineOffset1=IntegerPart[NoLatticePts/15]+1
];

Do[
ReadLine[Strm];
CurrLine++;
,{i,1,LineOffset1}];

(*We are now going to read the TB data onto a list:*)
DataLength=NoLatticePts*HamltnDim^2;

If[YesMessage==True,
Print["Wannier data begins at line ",CurrLine+1,"."];
];

LinesRead=0;
ReadAttempts=2;

Close[Strm];

TBData=Import[DataFile][[CurrLine+1;;]];

,
Print[Style["Invalid option specified for \"FileType\". Use \"CSV\" or \"Wannier90\".",Red]];
Return[{}]
]
];

MaxIndex=Max[Flatten[TBData[[All,SpaceDim+1;;SpaceDim+2]]]];

(*Print["Max orbital #: ",MaxIndex];*)

If[HamltnDim < MaxIndex,
HamltnDim = MaxIndex;
];

(*RowLength=Length[TBData[[1]]];*)

(*Print["Data shape: ",Dimensions[TBData]];*)

If[YesMessage==True,
Print["Tight-binding data file loaded."];
Print["Space dimension: ",SpaceDim];
Print["Number of orbitals: ",HamltnDim];
Print["Number of supercells: ",NoLatticePts];
Print["First row of data: ",TBData[[1]]];
Print["Last row of data: ",TBData[[-1]]];
Print["Number of entries: ",Length[TBData]];
Print["Loading Hamiltonian."];
];

kVec=Table[Symbol["k"<>ToString[i]],{i,1,SpaceDim}]; (*The symbolic k vector*)

If[Length[AtomLocn]==HamltnDim && AllTrue[AtomLocn,Length[#]==SpaceDim&] && AllTrue[Flatten[AtomLocn]//N,NumberQ[#]&],

If[OptionValue["SublatticePhase"]==True,
SublattPh=Table[(Cos[kVec . ((AtomLocn[[i]]-AtomLocn[[j]]) . BasisVectors)] + I Sin[kVec . ((AtomLocn[[i]]-AtomLocn[[j]]) . BasisVectors)]),{i,HamltnDim},{j,HamltnDim}];
,
SublattPh=1;
];

If[CheckSymmetrization[SymmList,BasisVectors,HamltnDim,SpaceDim,YesMessage,Resln]==True,
NeedsSymmetrization=True;
Normalisation=1/Length[SymmList];
If[YesMessage==True,
Print["Valid symmetrization scheme provided! We will symmetrize the loaded Hamiltonian."];
];

,
If[YesMessage==True,
Print[Style["No valid symmetrization scheme provided! We will not symmetrize the loaded Hamiltonian.",Red]];
];
NeedsSymmetrization=False;
]
,
If[YesMessage==True,
Print[Style["No valid symmetrization scheme provided! We will not symmetrize the loaded Hamiltonian.",Red]];
];
NeedsSymmetrization=False;
SublattPh=1;
];

(*The Hamiltonian matrix that we shall construct and return. And the symbolic k vector*)
Hamltn=ConstantArray[0,{HamltnDim,HamltnDim}];

LattVect=TBData[[1,1;;SpaceDim]]; (*The lattice vector that we iterate through in the data file. Initialised to the lattice vector of the first data entry.*)

TBMatEls={LattVect}; (*We shall contruct this list which contains the matrix elements in the data file,
indexed by the lattice vectors.*)
indx=1;

(*NoLines={Length[TBData]};*)

Do[
NewElmt=row[[SpaceDim+1;;SpaceDim+2]]; (*This is gives the ij indices corresponding the matrix element:
<i,R|H|j,0>*)
AppendTo[NewElmt,(row[[SpaceDim+3]]+I row[[SpaceDim+4]])];(*For each ij, include the matrix element drawn from the datafile so that the final form of
NewElmt is {i, j, Subscript[H, ij]}}*)

If[row[[1;;3]]==LattVect,
AppendTo[TBMatEls[[indx]],NewElmt]
,
LattVect=row[[1;;3]]; (*If a new lattice vector has been encountered, we have to create a entry.*)
indx++;
AppendTo[TBMatEls,LattVect];
AppendTo[TBMatEls[[indx]],NewElmt];
]
,{row,TBData}];

(*Sorting by the lattice vectors*)
TBMatEls=SortBy[TBMatEls,#[[1;;SpaceDim]]&]; (*Sorting the list by lattice vectors.*)

(*Sorting the data of each lattice point by the i,j index pair of the Hamiltonian matrix element*)
(*This is not needed for the tight binding loading function but needed for the interpolation routine.*)
(*Do[
TBMatEls[[i,SpaceDim+1;;Length[TBMatEls[[i]]]]]=SortBy[TBMatEls[[i,SpaceDim+1;;Length[TBMatEls[[i]]]]],#[[1;;2]]&]
(*Additionally, within each lattcice vector entry, we sort the elements by ij.*)
,{i,1,Length[TBMatEls]}];*)
(*Both these sorting procedures are needed to ensure that the binary search that we will perform below works properly.*)

LattPtsList=TBMatEls[[All,1;;SpaceDim]];

TBMat={}; (*This is a list of matrices with same dimension as the Hamiltonian, 
indexed by the lattice vector.*)

Do[
NewElmt=row[[1;;SpaceDim]]; 
AppendTo[NewElmt,ConstantArray[0,{HamltnDim,HamltnDim}]];

Do[

NewElmt[[SpaceDim+1,row[[i,1]],row[[i,2]]]]=row[[i,3]]; (*1 and 2 are ij index, 3 is the matrix element*)

,{i,SpaceDim+1,Length[row]}];

AppendTo[TBMat,NewElmt];

,{row,TBMatEls}];

If[NeedsSymmetrization,

TBMatExtras={};
(*We first iterate through the list of all the super cell lattice points and check if this list is closed under the given set of symmetries.*)
(*If not, we will insert a new element into TBMat and TBMatSymm at the appropriate location after some leg work*)
Do[
CurrLattPt=TBMat[[l,1;;SpaceDim]];

Do[
Do[
If[Chop[TBMat[[l,SpaceDim+1,i,j]],Resln]!=0,
Do[
NewLattPt=Floor[(Symm[[1]] . (CurrLattPt + AtomLocn[[j]]) + Symm[[2]])]-Floor[(Symm[[1]] . AtomLocn[[i]] + Symm[[2]])];

{LattPtLocn,Located}=BinSearch[LattPtsList,NewLattPt];
(*{LattPtLocn,Located}=BinSearch[TBMat[[All,1;;SpaceDim]],NewLattPt];*)

If[!Located,
(*This means that the supercell entry was not originally present in TBMat*)
(*We will now check if this new supercell has already been inserted into TBMatExtras in a previous iteration*)

If[Length[TBMatExtras]==0,
(*There were no supercells inserted until now so we should simply insert it*)
TBMatExtras={Append[NewLattPt,ConstantArray[0,{HamltnDim,HamltnDim}]]};
,
(*We will now check if this new supercell has already been inserted into TBMatExtras in a previous iteration*)
{LattPtLocn,Located}=BinSearch[TBMatExtras[[All,1;;SpaceDim]],NewLattPt];

If[!Located,
(*The supercell wasn't previously inserted into TBMatExtras so we should insert it now*)
TBMatExtras=Insert[TBMatExtras,Append[NewLattPt,ConstantArray[0,{HamltnDim,HamltnDim}]],LattPtLocn];
]
];
]

,{Symm,SymmList}]
]
,{j,1,HamltnDim}]
,{i,1,HamltnDim}]

,{l,Length[TBMat]}];

(*Print[TBMatExtras[[All,1;;SpaceDim]]];*)

(*Having created TBMatExtras with the extra supercells, we have to populate its matrix elements*)
(*We do this by iterating over the symmetries, then ij for each new supercell and add to H_ij appropriate linear combinations of 
the elements already existing in TBMat*)
(*This is necessary for otherwise, we may end up averaging some matrix elements with zero selectively in some directions 
during symmetrization*)
Do[
CurrLattPt=TBMatExtras[[l,1;;SpaceDim]];

Do[
Do[
NormalisationExtras=0;

Do[
NewLattPt=Floor[(Symm[[1]] . (CurrLattPt + AtomLocn[[j]]) + Symm[[2]])]-Floor[(Symm[[1]] . AtomLocn[[i]] + Symm[[2]])];

{LattPtLocn,Located}=BinSearch[LattPtsList,NewLattPt];
(*{LattPtLocn,Located}=BinSearch[TBMat[[All,1;;SpaceDim]],NewLattPt];*)

If[Located,
(*This means that the supercell entry was originally present in TBMat*)
(*We can use its elements to update H_ij in TBMatExtras*)
NormalisationExtras+=1;

TBMatExtras[[l,SpaceDim+1,i,j]]+=Simplify[((Symm[[3]] . Transpose[TBMat[[LattPtLocn,SpaceDim+1]]] . ConjugateTranspose[Symm[[3]]])[[j,i]])];
]

,{Symm,SymmList}];

If[NormalisationExtras!=0,
TBMatExtras[[l,SpaceDim+1,i,j]]=Simplify[1. TBMatExtras[[l,SpaceDim+1,i,j]] / NormalisationExtras];
]

,{j,1,HamltnDim}]
,{i,1,HamltnDim}];

,{l,Length[TBMatExtras]}];

(*Finally we can update TBMat using TBMatExtras*)
Do[
CurrLattPt=TBMatExtras[[l,1;;SpaceDim]];

{LattPtLocn,Located}=BinSearch[TBMat[[All,1;;SpaceDim]],CurrLattPt];

TBMat=Insert[TBMat,TBMatExtras[[l]],LattPtLocn];

,{l,Length[TBMatExtras]}];

(*We now update the list LattPtsList having inserted all the extra supercells into TBMat*)
LattPtsList=TBMat[[All,1;;SpaceDim]];

(*We now define the (to be) symmetrized version of TBMat*)
TBMatSymm=First@Last@Reap[
Do[
Sow[Append[LattPt,ConstantArray[0,{HamltnDim,HamltnDim}]]]
,{LattPt,LattPtsList}]];

Do[

Do[
CurrLattPt=TBMatSymm[[l,1;;SpaceDim]];

Do[
Do[
NewLattPt=Floor[(Symm[[1]] . (CurrLattPt + AtomLocn[[j]]) + Symm[[2]])]-Floor[(Symm[[1]] . AtomLocn[[i]] + Symm[[2]])];

{LattPtLocn,Located}=BinSearch[LattPtsList,NewLattPt];

If[Located,
TBMatSymm[[l,SpaceDim+1,i,j]]+=Simplify[((Symm[[3]] . Transpose[TBMat[[LattPtLocn,SpaceDim+1]]] . ConjugateTranspose[Symm[[3]]])[[j,i]])Normalisation];
(*Instead of the convoluted expression above, we might like to use the following:*)
(*TBMatSymm[[l,SpaceDim+1,i,j]]+=Simplify[((Symm[[3]] . TBMat[[LattPtLocn,SpaceDim+1]] . ConjugateTranspose[Symm[[3]]])[[i,j]])Normalisation];*)
(*But we don't use this since Mathematica defines matrix multiplication in a peculiar way.*)

TBMatSymm[[l,SpaceDim+1,i,j]]=Simplify[TBMatSymm[[l,SpaceDim+1,i,j]]];

(*
,
If[TBMat[[BinSearch[LattPtsList,CurrLattPt][[1]],SpaceDim+1,i,j]]!=0,
Print[NewLattPt,Style[" Not Located",Red]," for CurrPt = ",CurrLattPt," and (i,j) = (",i,",",j,") with Hij = ",TBMat[[BinSearch[LattPtsList,CurrLattPt][[1]],SpaceDim+1,i,j]]];
]*)
]
,{j,1,HamltnDim}]
,{i,1,HamltnDim}]


,{l,1,Length[LattPtsList]}]

,{Symm,SymmList}]
,
TBMatSymm=TBMat;
];

If[YesMessage==True && NeedsSymmetrization==True,
Print["Number of supercells after symmetrization: ",Length[TBMatSymm]];
];

(*Finally, we shall construct the Hamiltonian*)
Do[
Hamltn+=Simplify[0.5(Chop[TBRow[[SpaceDim+1]],Resln] SublattPh (Cos[kVec . (TBRow[[1;;SpaceDim]] . BasisVectors)] + I Sin[kVec . (TBRow[[1;;SpaceDim]] . BasisVectors)])+ ComplexExpand[ConjugateTranspose[Chop[TBRow[[SpaceDim+1]],Resln] SublattPh]](Cos[kVec . (TBRow[[1;;SpaceDim]] . BasisVectors)] - I Sin[kVec . (TBRow[[1;;SpaceDim]] . BasisVectors)]))];
,{TBRow,TBMatSymm}];

{Hamltn,TBMatSymm}
]


(* ::Section::Closed:: *)
(*Loading and interpolating through a series of tight binding models stored in Wannier90 format, with symmetrization*)


Options[LoadInterpolateTBMs]={"FileType"->"Wannier90","LatticeVectors"->-1,"Dimension"->-1,"SpaceDimension"->3,"Messages"->True,"OrbitalCentres"->-1,"Symmetries"->-1,"Resolution"->10^-10,"WannierFileOffset"->4,"SublatticePhase"->False,"ParameterSymbol"->"\[Delta]t","PolynomialDegree"->-1}

(*This funtion loads a tight binding model, stored in a csv or space sparated data file into a Hamiltonian.*)
LoadInterpolateTBMs[DataFiles_,ParamList_,OptionsPattern[]]:=Module[
{Strm,StrmList,Monoms,SkipLines,ParamSymbol,TBMatInterpolated,IndxLocn,InterpolationDegree,HopElmt,HopFitData,NoParams,SublattPh,CurrLine,TempLine,NoLatticePts,ReadAttempts,LinesRead,DataLength,TBMatExtras,NormalisationExtras,LineOffset1,MaxLineOffset,BasisVectors,LattPtsList,LattPtLocn,NewLattPt,CurrLattPt,Located,TBMat,TBMatEls,LattVect,indx,NoLines,NewElmt,NeedsSymmetrization,SymmList,AtomLocn,Normalisation,TBMatInterpolatedSymm,TBData,Hamltn,kVec,MaxIndex,RowLength,HamltnDimList,HamltnDim,SpaceDim,YesMessage,Resln},

If[Length[DataFiles]==Length[ParamList]&&AllTrue[DataFiles,StringQ[#]&]&&AllTrue[ParamList,NumberQ[#]&],

SpaceDim=OptionValue["SpaceDimension"];
(*YesMessage=OptionValue["Messages"];*)
(*HamltnDim=OptionValue["Dimension"];*)
If[BooleanQ[OptionValue["Messages"]],
YesMessage=OptionValue["Messages"];
,
YesMessage=True;
];

AtomLocn=OptionValue["OrbitalCentres"];
SymmList=OptionValue["Symmetries"];
Resln=OptionValue["Resolution"];
MaxLineOffset=OptionValue["WannierFileOffset"];
ParamSymbol=OptionValue["ParameterSymbol"];

NoParams=Length[DataFiles];
HamltnDimList=ConstantArray[0,NoParams];
NoLatticePts=ConstantArray[0,NoParams];
SkipLines=ConstantArray[0,NoParams];
StrmList=ConstantArray[0,NoParams];

If[IntegerQ[OptionValue["PolynomialDegree"]]&&OptionValue["PolynomialDegree"]>0,
InterpolationDegree=OptionValue["PolynomialDegree"];
,
InterpolationDegree=Length[ParamList]-1
(*The degree to be used for the polynomial interpolation.*)
];

,
Print[Style["Invalid input(s)! The first input variable should be a list of strings containing the Wannier90 datafile adresses and the second variable should be a list of numbers, containing the values of the corresponding parameters.",Red]];
Return[{}]
];

(*Print[OptionValue["LatticeVectors"]];*)

If[OptionValue["LatticeVectors"]==-1||Length[OptionValue["LatticeVectors"]]!=SpaceDim||!AllTrue[Flatten[OptionValue["LatticeVectors"]]//N,NumberQ[#]&]||!AllTrue[OptionValue["LatticeVectors"],Length[#]==SpaceDim&],
(*If no basis vectors have been provided or if invalid ones are provided, we assume that it is square or cubic lattice*)
If[YesMessage==True,
Print["Invalid or no lattice vectors provided! Using unit vectors."];
];

BasisVectors=IdentityMatrix[SpaceDim];
,
If[YesMessage==True,
Print["Valid lattice vectors provided: ",BasisVectors];
];

BasisVectors=OptionValue["LatticeVectors"];
];

(*Print[BasisVectors];*)

Do[
StrmList[[j]]=OpenRead[DataFiles[[j]]];

(*Skipping the first few lines that just contain the record of when the file was written.*)
(*The variable SkipLines will contain the number of such lines in each file.*)
Do[
TempLine=ToExpression[StringSplit[ReadLine[StrmList[[j]]]][[1]]];

SkipLines[[j]]++;

If[IntegerQ[TempLine],
HamltnDimList[[j]]=TempLine;
Break[];
]
,{i,1,MaxLineOffset}];

(*Reading the number of lattice points over which the hopping terms are specified.*)
NoLatticePts[[j]]=ToExpression[ReadLine[StrmList[[j]]]];
SkipLines[[j]]++;

(*Calculating the number lines which have the form 1 1 1 1 1 1 ...1 1 1.*)
If[IntegerQ[NoLatticePts[[j]]/15],
SkipLines[[j]]+=IntegerPart[NoLatticePts[[j]]/15];
,
SkipLines[[j]]+=IntegerPart[NoLatticePts[[j]]/15]+1
];

Close[StrmList[[j]]];

,{j,1,NoParams}];

HamltnDim=Max[HamltnDimList];

(*Although the code is designed to interpolate between data-files/models with different number of orbitals and supercells, we disallow this feature in the current version since it has not been extensively benchmarked.*)
If[!AllTrue[HamltnDimList,#==HamltnDim&]||!AllTrue[NoLatticePts,#==NoLatticePts[[1]]&],
(*If the number of orbitals and supercells are not verified above to be the same for all the Wannier90 models, the data files are not compatible and the hopping matrix elements cannot be interpolated with a function.*)
Print[Style["Incompatible data files! All the Wannier90 data files should have the same number of orbitals and super cells.",Red]];
Return[{}]
];

DataLength=NoLatticePts[[1]]*HamltnDim^2;

(*Loading the first file. This will dictate the basic structure of the data.*)
TBData=Import[DataFiles[[1]]][[SkipLines[[1]]+1;;]];
LattVect=TBData[[1,1;;SpaceDim]]; (*The lattice vector that we iterate through in the data file.*)
TBMatEls={LattVect}; (*We shall contruct this list which contains the matrix elements for each of the data files,
indexed by the lattice vectors.*)
indx=1;
NoLines={Length[TBData]};

Do[
NewElmt=row[[SpaceDim+1;;SpaceDim+2]]; (*This is gives the ij indices corresponding the matrix element:
<i,R|H|j,0>*)
AppendTo[NewElmt,Join[{row[[SpaceDim+3]]+I row[[SpaceDim+4]]},ConstantArray[0,NoParams-1]]];(*For each ij, we 
include a list of matrix elements drawn from the datafiles that we shall interpolate through. The final form of
NewElmt is {i, j, {\!\(
\*SubsuperscriptBox[\(H\), \(ij\), \((1)\)], \ 
\*SubsuperscriptBox[\(H\), \(ij\), \((2)\)], \ \(\(...\)\(..\)\)\)}}*)

If[row[[1;;3]]==LattVect,
AppendTo[TBMatEls[[indx]],NewElmt]
,
LattVect=row[[1;;3]]; (*If a new lattice vector has been encountered, we have to create a entry.*)
indx++;
AppendTo[TBMatEls,LattVect];
AppendTo[TBMatEls[[indx]],NewElmt];
]
,{row,TBData}];

(*Sorting by the lattice vectors*)
TBMatEls=SortBy[TBMatEls,#[[1;;SpaceDim]]&]; (*Sorting the list by lattice vector.*)

(*Sorting the data of each lattice point by the i,j index pair of the Hamiltonian matrix element*)
Do[
TBMatEls[[i,SpaceDim+1;;Length[TBMatEls[[i]]]]]=SortBy[TBMatEls[[i,SpaceDim+1;;Length[TBMatEls[[i]]]]],#[[1;;2]]&]
(*Additionally, within each lattcice vector entry, we sort the elements by ij.*)
,{i,1,Length[TBMatEls]}];

LattPtsList=TBMatEls[[All,1;;SpaceDim]];

(*Loading the rest of the data files. Each new entry has to be checked for existence.*)
Do[
TBData=Import[DataFiles[[i]]][[SkipLines[[i]]+1;;]];

AppendTo[NoLines,Length[TBData]];

Do[
{LattPtLocn,Located}=BinSearch[LattPtsList,row[[1;;SpaceDim]]];

If[!Located,
(*The lattice point was not located. We need to insert a new entry.*)
TBMatEls=Insert[TBMatEls,Append[row[[1;;SpaceDim]](*the lattice vector*),Append[row[[SpaceDim+1;;SpaceDim+2]](*the i-j index*),Join[ConstantArray[0,i-1],ConstantArray[row[[SpaceDim+3]]+I row[[SpaceDim+4]],1],ConstantArray[0,NoParams-i]](*the data array*)]],LattPtLocn]
,
(*The lattice point has been located. We need to check if the entry for the ij index is already there.*)

{IndxLocn,Located}=BinSearch[TBMatEls[[LattPtLocn,SpaceDim+1;;,1;;2]],row[[SpaceDim+1;;SpaceDim+2]]];

If[!Located,
(*Although the lattice point was located, the particular ij data was not. We therefore need to insert a new element for the ij index.*)
TBMatEls=Insert[TBMatEls,Append[row[[SpaceDim+1;;SpaceDim+2]](*the i-j index*),Join[ConstantArray[0,i-1],ConstantArray[row[[SpaceDim+3]]+I row[[SpaceDim+4]],1],ConstantArray[0,NoParams-i]]],{LattPtLocn,SpaceDim+IndxLocn}];
,
(*The index also has been located. So no insertion is needed. We merely update the data at the appropriate location.*)
TBMatEls[[LattPtLocn,SpaceDim+IndxLocn,3,i]]=row[[SpaceDim+3]]+I row[[SpaceDim+4]];
]
]
,{row,TBData}];
,{i,2,NoParams}];


If[YesMessage==True,
Print["Tight-binding data files loaded."];
Print["Space dimension: ",SpaceDim];
Print["Number of orbitals: ",HamltnDim];
Print["Number of supercells: ",NoLatticePts[[1]]];
Print["Number of models/data-files to interpolate: ",NoParams];
Print["Number of lines in the data files: ",(Table[SkipLines[[i]]+NoLines[[i]],{i,1,NoParams}])];
(*Print["First row: ",TBMatEls[[1]]];Print["Last row: ",TBMatEls[[-1]]];*)
Print["Interpolating between the Hamiltonians."];
];

kVec=Table[Symbol["k"<>ToString[i]],{i,1,SpaceDim}]; (*The symbolic k vector*)
(*We need to generate the set of monomials to be used for the fit*)

Monoms=First@Last@Reap[
Do[
Sow[(Symbol[ParamSymbol])^i]
,{i,0,InterpolationDegree}]
];

(*We check if a valid symmetrization scheme is provided and if the sublattice phase is to be included*)
If[Length[AtomLocn]==HamltnDim && AllTrue[AtomLocn,Length[#]==SpaceDim&] && AllTrue[Flatten[AtomLocn]//N,NumberQ[#]&],

If[OptionValue["SublatticePhase"]==True,
SublattPh=Table[(Cos[kVec . ((AtomLocn[[i]]-AtomLocn[[j]]) . BasisVectors)] + I Sin[kVec . ((AtomLocn[[i]]-AtomLocn[[j]]) . BasisVectors)]),{i,HamltnDim},{j,HamltnDim}];
,
SublattPh=1;
];

If[CheckSymmetrization[SymmList,BasisVectors,HamltnDim,SpaceDim,YesMessage,Resln]==True,
NeedsSymmetrization=True;
Normalisation=1/Length[SymmList];
If[YesMessage==True,
Print["Valid symmetrization scheme provided! We will symmetrize the loaded Hamiltonian."];
];

,
If[YesMessage==True,
Print[Style["No valid symmetrization scheme provided! We will not symmetrize the loaded Hamiltonian.",Red]];
];
NeedsSymmetrization=False;
]
,
If[YesMessage==True,
Print[Style["No valid symmetrization scheme provided! We will not symmetrize the loaded Hamiltonian.",Red]];
];
NeedsSymmetrization=False;
SublattPh=1;
];

(*The Hamiltonian matrix that we shall construct and return. And the symbolic k vector*)
Hamltn=ConstantArray[0,{HamltnDim,HamltnDim}];

TBMatInterpolated={}; (*This is a list of matrices with same dimension as the Hamiltonian, 
indexed by the lattice vector.*)

Do[
NewElmt=row[[1;;SpaceDim]];
AppendTo[NewElmt,ConstantArray[0,{HamltnDim,HamltnDim}]];

Do[

(*Data set for polynomial fit*)
HopFitData=First@Last@Reap[Do[
Sow[{ParamList[[j]],Chop[row[[i,3(*1 and 2 are ij index*),j]],Resln]}];
,{j,1,NoParams}]];

(*Polynomial fit for the hopping element*)
HopElmt=Fit[HopFitData,Monoms,Symbol[ParamSymbol]];

NewElmt[[SpaceDim+1,row[[i,1]],row[[i,2]]]]=HopElmt;
(*
HopElmtConj=ComplexExpand[Conjugate[HopElmt]];
(*=Fit[Conjugate[HopFitData],Monoms,Symbol[ParamSymbol]];*)

Hamltn[[row[[i,1]],row[[i,2]]]]+=0.5 HopElmt(Cos[kVec . row[[1;;SpaceDim]]]-I Sin[kVec . row[[1;;SpaceDim]]]);
Hamltn[[row[[i,2]],row[[i,1]]]]+=0.5 HopElmtConj(Cos[kVec . row[[1;;SpaceDim]]]+I Sin[kVec . row[[1;;SpaceDim]]]);
*)
,{i,SpaceDim+1,Length[row]}];

AppendTo[TBMatInterpolated,NewElmt];

,{row,TBMatEls}];

(*Now we symmetrize TBMatInterpolated.*)
(*The list SymmList contains all the necessary symmetries of the space group. It is organised as follows:
Each element Symm of SymmList corresponds to a space group symmetry g. 
Symm[[1]] denotes the matrix for the orthogonal transformation part of g.
Symm[[2]] deontes the translation part of g.
Symm[[3]] deontes the representation of orthogonal transformation part on the Hilbert space, which is a 
unitary matrix of dimension HamltnDim.
As we iterate through lattice sites and i, j at each lattice site:
LatticeSite -> Floor[Symm[[1]].(LatticeSite + AtomLocn[[i]])+Symm[[2]]]-Floor[Symm[[1]].AtomLocn[[j]]+Symm[[2]]],
TBMatInterpolated[[Indx[LatticeSite],SpaceDim+1]]-> ConjugateTranspose[Symm[[3]]].TBMatInterpolated[[Indx[NewLatticeSite],SpaceDim+1]].Symm[[3]]
*)

If[NeedsSymmetrization,

TBMatExtras={};
(*We first iterate through the list of all the super cell lattice points and check if this list is closed under the given set of symmetries.*)
(*If not, we will insert a new element into TBMat and TBMatSymm at the appropriate location after some leg work*)
Do[
CurrLattPt=TBMatInterpolated[[l,1;;SpaceDim]];

Do[
Do[
If[!PossibleZeroQ[Chop[TBMatInterpolated[[l,SpaceDim+1,i,j]],Resln]],
Do[
NewLattPt=Floor[(Symm[[1]] . (CurrLattPt + AtomLocn[[j]]) + Symm[[2]])]-Floor[(Symm[[1]] . AtomLocn[[i]] + Symm[[2]])];

{LattPtLocn,Located}=BinSearch[LattPtsList,NewLattPt];
(*{LattPtLocn,Located}=BinSearch[TBMat[[All,1;;SpaceDim]],NewLattPt];*)

If[!Located,
(*This means that the supercell entry was not originally present in TBMat*)
(*We will now check if this new supercell has already been inserted into TBMatExtras in a previous iteration*)

If[Length[TBMatExtras]==0,
(*There were no supercells inserted until now so we should simply insert it*)
TBMatExtras={Append[NewLattPt,ConstantArray[0,{HamltnDim,HamltnDim}]]};
,
(*We will now check if this new supercell has already been inserted into TBMatExtras in a previous iteration*)
{LattPtLocn,Located}=BinSearch[TBMatExtras[[All,1;;SpaceDim]],NewLattPt];

If[!Located,
(*The supercell wasn't previously inserted into TBMatExtras so we should insert it now*)
TBMatExtras=Insert[TBMatExtras,Append[NewLattPt,ConstantArray[0,{HamltnDim,HamltnDim}]],LattPtLocn];
]
];
]

,{Symm,SymmList}]
]
,{j,1,HamltnDim}]
,{i,1,HamltnDim}]

,{l,Length[TBMatInterpolated]}];

(*Print[TBMatExtras[[All,1;;SpaceDim]]];*)

(*Having created TBMatExtras with the extra supercells, we have to populate its matrix elements*)
(*We do this by iterating over the symmetries, then ij for each new supercell and add to H_ij appropriate linear combinations of 
the elements already existing in TBMat*)
(*This is necessary for otherwise, we may end up averaging some matrix elements with zero selectively in some directions 
during symmetrization*)
(*
Do[
CurrLattPt=TBMatExtras[[l,1;;SpaceDim]];

Do[
Do[
NormalisationExtras=0;

Do[
NewLattPt=Floor[(Symm[[1]] . (CurrLattPt + AtomLocn[[j]]) + Symm[[2]])]-Floor[(Symm[[1]] . AtomLocn[[i]] + Symm[[2]])];

{LattPtLocn,Located}=BinSearch[LattPtsList,NewLattPt];
(*{LattPtLocn,Located}=BinSearch[TBMat[[All,1;;SpaceDim]],NewLattPt];*)

If[Located,
(*This means that the supercell entry was originally present in TBMat*)
(*We can use its elements to update H_ij in TBMatExtras*)
NormalisationExtras+=1;

TBMatExtras[[l,SpaceDim+1,i,j]]+=Simplify[((Symm[[3]] . Transpose[TBMatInterpolated[[LattPtLocn,SpaceDim+1]]] . ConjugateTranspose[Symm[[3]]])[[j,i]])];
]

,{Symm,SymmList}];

If[NormalisationExtras!=0,
TBMatExtras[[l,SpaceDim+1,i,j]]=Simplify[1. TBMatExtras[[l,SpaceDim+1,i,j]] / NormalisationExtras];
]

,{j,1,HamltnDim}]
,{i,1,HamltnDim}];

,{l,Length[TBMatExtras]}];
*)

(*Finally we can update TBMat using TBMatExtras*)
Do[
CurrLattPt=TBMatExtras[[l,1;;SpaceDim]];

{LattPtLocn,Located}=BinSearch[TBMatInterpolated[[All,1;;SpaceDim]],CurrLattPt];

TBMatInterpolated=Insert[TBMatInterpolated,TBMatExtras[[l]],LattPtLocn];

,{l,Length[TBMatExtras]}];

(*We now update the list LattPtsList having inserted all the extra supercells into TBMat*)
LattPtsList=TBMatInterpolated[[All,1;;SpaceDim]];


TBMatInterpolatedSymm=First@Last@Reap[
Do[
Sow[Append[LattPt,ConstantArray[0,{HamltnDim,HamltnDim}]]]
,{LattPt,LattPtsList}]];

Do[

Do[
CurrLattPt=TBMatInterpolatedSymm[[l,1;;SpaceDim]];

Do[
Do[
NewLattPt=Floor[Symm[[1]] . (CurrLattPt + AtomLocn[[i]])+Symm[[2]]]-Floor[Symm[[1]] . AtomLocn[[j]]+Symm[[2]]];

{LattPtLocn,Located}=BinSearch[LattPtsList,NewLattPt];

If[Located,
TBMatInterpolatedSymm[[l,SpaceDim+1,i,j]]+=Simplify[((Symm[[3]] . Transpose[TBMatInterpolated[[LattPtLocn,SpaceDim+1]]] . ConjugateTranspose[Symm[[3]]])[[j,i]])Normalisation];
TBMatInterpolatedSymm[[l,SpaceDim+1,i,j]]=Simplify[TBMatInterpolatedSymm[[l,SpaceDim+1,i,j]]]
(*
,
If[TBMat[[BinSearch[LattPtsList,CurrLattPt][[1]],SpaceDim+1,i,j]]!=0,
Print[NewLattPt,Style[" Not Located",Red]," for CurrPt = ",CurrLattPt," and (i,j) = (",i,",",j,") with Hij = ",TBMat[[BinSearch[LattPtsList,CurrLattPt][[1]],SpaceDim+1,i,j]]];
]
*)
]
,{j,1,HamltnDim}]
,{i,1,HamltnDim}]


,{l,1,Length[LattPtsList]}]

,{Symm,SymmList}]
,

TBMatInterpolatedSymm=TBMatInterpolated;
];

If[YesMessage==True && NeedsSymmetrization==True,
Print["Number of supercells after symmetrization: ",Length[TBMatInterpolatedSymm]];
];

(*Finally, we shall construct the Hamiltonian*)
Do[
(*Hamltn+=Simplify[0.5(Chop[TBRow[[SpaceDim+1]],Resln] SublattPh (Cos[kVec . (TBRow[[1;;SpaceDim]] . BasisVectors)] + I Sin[kVec . (TBRow[[1;;SpaceDim]] . BasisVectors)])+ ComplexExpand[ConjugateTranspose[Chop[TBRow[[SpaceDim+1]],Resln] SublattPh]](Cos[kVec . (TBRow[[1;;SpaceDim]] . BasisVectors)] - I Sin[kVec . (TBRow[[1;;SpaceDim]] . BasisVectors)]))];
*)
(*TBRow[[1;;SpaceDim]]*)
Hamltn+=Simplify[0.5(TBRow[[SpaceDim+1]](Cos[kVec . (TBRow[[1;;SpaceDim]] . BasisVectors)] + I Sin[kVec . (TBRow[[1;;SpaceDim]] . BasisVectors)])+ComplexExpand[ConjugateTranspose[TBRow[[SpaceDim+1]]]](Cos[kVec . (TBRow[[1;;SpaceDim]] . BasisVectors)] - I Sin[kVec . (TBRow[[1;;SpaceDim]] . BasisVectors)]))];
,{TBRow,TBMatInterpolatedSymm}];

{Hamltn,TBMatInterpolatedSymm,ParamSymbol,ParamList}

]


(* ::Section::Closed:: *)
(*Writing a TBM Hamiltonian to a CPP file*)


Options[WriteHamiltonianCPP]={"FilePath"->-1,"Messages"->True,"ParameterDimension"->0}

WriteHamiltonianCPP[Hamltn_,SpaceDim_,OptionsPattern[]]:=
Module[
{FileNamePath,ParamDim,YesMessage,strm,strm2,kVec,tVec,kVecFull,HamltnDim},

If[BooleanQ[OptionValue["Messages"]],
YesMessage=OptionValue["Messages"];
,
YesMessage=True;
];

(*We begin by checking if a valid filename and/or a valid directory have been provided*)
(*If not, we use a default name*)
If[!StringQ[OptionValue["FilePath"]],
FileNamePath=NotebookDirectory[]<>"cpp/Hamltn.cpp";

If[YesMessage==True,
Print[Style["Invalid or no file path provided! Reverting to the default name Hamltn.cpp. File shall be stored in the directory cpp located in the local directory of the notebook.",Red]];
]

,

If[DirectoryQ[StringDrop[OptionValue["FilePath"],{StringPosition[OptionValue["FilePath"],"/"][[-1,1]]+1,-1}]] && StringDrop[OptionValue["FilePath"],StringPosition[OptionValue["FilePath"],"."][[-1,1]]]=="cpp",
FileNamePath=OptionValue["FilePath"];

,

FileNamePath=NotebookDirectory[]<>"cpp/Hamltn.cpp";

If[!DirectoryQ[NotebookDirectory[]<>"cpp"],
CreateDirectory[NotebookDirectory[]<>"cpp"]
];

If[YesMessage==True,
Print[Style["Invalid or no file path or extension provided! Reverting to the default name Hamltn.cpp. File shall be stored in the directory cpp located in the local directory of the notebook.",Red]];
]
]
];

(*We now check if a valid space dimension and parameter dimension have been provided*)
If[IntegerQ[SpaceDim] && SpaceDim>0,

If[IntegerQ[OptionValue["ParameterDimension"]] && 0 <= OptionValue["ParameterDimension"] < SpaceDim,

ParamDim=OptionValue["ParameterDimension"];

,

ParamDim=0;

If[YesMessage==True,
Print[Style["Invalid parameter dimension! It has to be a non-negative integer less than the space dimension. Setting it to the default value of zero.",Red]]
]
];

If[AllTrue[Flatten[Hamltn[ConstantArray[0,SpaceDim+ParamDim]]//N],NumberQ[#]&] && SquareMatrixQ[Hamltn[ConstantArray[0,SpaceDim+ParamDim]]],
(*Everything is fine... so far*)
strm=OpenWrite[FileNamePath];

kVec=Table[Symbol["k"<>ToString[i]],{i,1,SpaceDim}];
tVec=Table[Symbol["t"<>ToString[i]],{i,1,ParamDim}];
kVecFull=Join[kVec,tVec];

HamltnDim=Length[Hamltn[ConstantArray[0,SpaceDim+ParamDim]]];

(*We can now start writing the Hamiltonian in the file*)
WriteString[strm,"const int hmlt_dim = "<>ToString[HamltnDim]<>";\n"];
WriteString[strm,"const int param_dim = "<>ToString[ParamDim]<>";\n"];
WriteString[strm,"const int vec_dim = "<>ToString[SpaceDim]<>";\n\n"];

WriteString[strm,"void hmlt_def(complex<double> **hmlt"];

Do[
WriteString[strm,", double "<>ToString[ki]];
,{ki,kVecFull}];

WriteString[strm,")\n"];

WriteString[strm,"{\n"];

Do[
Do[
WriteString[strm,"\t"<>"hmlt["<>ToString[i-1]<>"]["<>ToString[j-1]<>"] = "];

WriteString[strm,StringDelete[ToString[CForm[Simplify[Hamltn[kVecFull][[i,j]]]//N]],"\n"|"\r"]];

WriteString[strm,";\n"]
,{j,i,HamltnDim}]
,{i,1,HamltnDim}];

WriteString[strm,"\n\n"];

Do[
Do[
WriteString[strm,"\t"<>"hmlt["<>ToString[i-1]<>"]["<>ToString[j-1]<>"] = conj("<>"hmlt["<>ToString[j-1]<>"]["<>ToString[i-1]<>"]);\n"];
,{j,1,i-1}];
,{i,2,HamltnDim}];

WriteString[strm,"}"];
Close[strm];

(*We now create another file containing the necessary functions for the C++ code to be able to compile properly.*)
strm2=OpenWrite[StringDrop[FileNamePath,{StringPosition[FileNamePath,"/"][[-1,1]]+1,-1}]<>"ancillaries.cpp"];

WriteString[strm2,"const complex<double> I(0.0,1.0);"<>"\n"<>"\n"];
WriteString[strm2,"// We begin by (re)defining some functions to accommodate Mathematica's sort of unusual CForm convention.\n"];

WriteString[strm2,"inline complex<double> Complex(double a, double b)\n"];
WriteString[strm2,"{\n"];
WriteString[strm2,"\t"<>"return (1. * a + 1. * b * I);\n"];
WriteString[strm2,"}\n"<>"\n"];

WriteString[strm2,"inline double Cos(double a)\n"];
WriteString[strm2,"{\n"];
WriteString[strm2,"\t"<>"return cos(a);\n"];
WriteString[strm2,"}\n"<>"\n"];

WriteString[strm2,"inline double Power(double a, double b)\n"];
WriteString[strm2,"{\n"];
WriteString[strm2,"\t"<>"return pow(a, b);\n"];
WriteString[strm2,"}\n"<>"\n"];

Close[strm2];

If[YesMessage==True,
Print["Hamiltonian successfully written to ",StringDrop[FileNamePath,StringPosition[FileNamePath,"/"][[-1,1]]],". Please include both ",StringDrop[FileNamePath,StringPosition[FileNamePath,"/"][[-1,1]]]," and ancillaries.cpp in your code."];
]
,

If[YesMessage==True,
Print[Style["Invalid Hamiltonian provided!",Red]]
]
]

,
If[YesMessage==True,
Print[Style["Space dimension has to be a positive integer!",Red]]
]
]
]


(* ::Section::Closed:: *)
(*Writing a plain or interpolated tight binding data to a file in Wannier90 format*)


(* ::Code::Initialization::Plain:: *)
Options[WriteTBM]={"ParameterSymbol"->-1,"ParameterValue"->0,"SpinOrbitMatrix"->-1,"BandRenormalization"->1,"FermiLevel"->0.,"Resolution"->10^-8,"Messages"->True,"CheckTightBindingData"->True}

WriteTBM[WannierFile_,TBData_,OptionsPattern[]]:=Module[
{file,strm,SkipLines,OneString,ParamSymbol,ParamVal,ValidSym,ResidueOnes,HoppMat,HoppMat2,TBRow,RowIndx,CurrRow,HopFitData,HoppFitElmt,NoDigits,NoLatticePts,CurrLattPt,ZRenorm,EFermi,SOC,ValidSOC,Resln,SpaceDim,HamltnDim,AllFine,YesMessage},

AllFine=True;

If[BooleanQ[OptionValue["Messages"]],
YesMessage=OptionValue["Messages"];
,
YesMessage=True;
];

If[StringQ[OptionValue["ParameterSymbol"] && !StringContainsQ[OptionValue["ParameterSymbol"]," "] && LetterQ[StringPart[OptionValue["ParameterSymbol"],1]]],
ParamSymbol=OptionValue["ParameterSymbol"];
ValidSym=True;
,
ValidSym=False;
];

If[NumberQ[OptionValue["ParameterValue"]] && OptionValue["ParameterValue"]\[Element]Reals,
ParamVal=OptionValue["ParameterValue"];
,
ParamVal=0;
];

(*To save time the user may have disabled the option for checking the format tight binding data, which may be very large.*)
If[!BooleanQ[OptionValue["CheckTightBindingData"]] || (BooleanQ[OptionValue["CheckTightBindingData"]] && OptionValue["CheckTightBindingData"]==True),
If[!ListQ[TBData] || !AllTrue[TBData,Length[#]==Length[TBData[[1]]]&] || Length[TBData[[1]]]<=1 || !AllTrue[TBData,Dimensions[#[[-1]]]=={Length[TBData[[1,-1]]],Length[TBData[[1,-1]]]}&] || !AllTrue[Flatten[TBData[[All,1;;-2]]],IntegerQ[#]&] || (ValidSym && !AllTrue[Flatten[ReplaceAll[TBData[[All,-1]],{Symbol[ParamSymbol]->ParamVal}]//Simplify],NumberQ[#]&]),
If[YesMessage==True,
Print[Style["Invalid tight binding data provided!",Red]];
];
AllFine=False;
Return[Null]
]
];

SpaceDim=Length[TBData[[1]]]-1;
HamltnDim=Length[TBData[[1,-1]]];

If[!MatrixQ[OptionValue["SpinOrbitMatrix"]] || Dimensions[OptionValue["SpinOrbitMatrix"]]!={2HamltnDim,2HamltnDim} || !AllTrue[Flatten[OptionValue["SpinOrbitMatrix"]]//N,NumberQ[#]&] || !HermitianMatrixQ[OptionValue["SpinOrbitMatrix"]],
ValidSOC=False;
,
ValidSOC=True;
SOC=OptionValue["SpinOrbitMatrix"];
];

If[NumberQ[OptionValue["FermiLevel"]] && OptionValue["FermiLevel"]\[Element]Reals,
EFermi=OptionValue["FermiLevel"];
,
EFermi=0.;

If[YesMessage,
Print["Invalid Fermi level provided. Setting it to zero."];
]
];

If[NumberQ[OptionValue["Resolution"]] && 0 < OptionValue["Resolution"] < 1,
Resln=OptionValue["Resolution"];
,
Resln=10^-8;

If[YesMessage,
Print["Invalid resolution provided. Reverting to the default value of \!\(\*SuperscriptBox[\(10\), \(-8\)]\)."]
]
];

If[NumberQ[OptionValue["BandRenormalization"]] && Chop[OptionValue["BandRenormalization"],Resln] > 0,
ZRenorm=OptionValue["BandRenormalization"];
,
ZRenorm=1.;

If[YesMessage,
Print["Invalid overall band renormalization provided. Setting it to 1."];
]
];

(*We begin by checking if a valid filename and/or a valid directory have been provided*)
(*If not, we use a default name*)
If[!StringQ[WannierFile],
file=NotebookDirectory[]<>"Wannier90/wannier90_hr.dat";

If[!DirectoryQ[NotebookDirectory[]<>"Wannier90"],
CreateDirectory[NotebookDirectory[]<>"Wannier90"]
];

If[YesMessage==True,
Print[Style["Invalid file path provided! Reverting to the default name wannier90_hr.dat. The file will be stored in the directory Wannier90 located in the local directory of the notebook.",Red]];
]
,

If[DirectoryQ[StringDrop[WannierFile,{StringPosition[WannierFile,"/"][[-1,1]]+1,-1}]] && StringDrop[WannierFile,StringPosition[WannierFile,"."][[-1,1]]]=="dat",
file=WannierFile;
,

file=NotebookDirectory[]<>"Wannier90/wannier90_hr.dat";

If[!DirectoryQ[NotebookDirectory[]<>"Wannier90"],
CreateDirectory[NotebookDirectory[]<>"Wannier90"]
];

If[YesMessage==True,
Print[Style["Invalid file path or extension provided! Reverting to the default name wannier90_hr.dat. The file will be stored in the directory Wannier90 located in the local directory of the notebook.",Red]];
]
]
];


If[AllFine,

NoDigits=-IntegerPart[Log10[Resln]];

NoLatticePts=Length[TBData];

(*file=NotebookDirectory[]<>WannierFile;*)

strm=OpenWrite[file];

(*We begin by writing the header, which contains information regarding the date and time of writing, number of orbitals and number of lattice points.*)WriteString[strm," written on "<>DateString[{"Day","MonthName","Year"}]<>" at "<>DateString[{"Hour",":","Minute",":","Second"}]<>"\r\n"];

If[ValidSOC,

WriteString[strm,StringPadLeft[ToString[2HamltnDim],12]<>"\r\n"]
,

WriteString[strm,StringPadLeft[ToString[HamltnDim],12]<>"\r\n"]
];

WriteString[strm,"          "<>ToString[NoLatticePts]<>"\r\n"];

(*We now proceed to write the strings containing ones*)If[IntegerQ[NoLatticePts/15],
SkipLines=IntegerPart[NoLatticePts/15];

OneString="";
Do[
OneString=StringJoin[OneString,"    1"]
,{j,1,15}];

Do[
WriteString[strm,OneString<>"\r\n"];
,{i,1,SkipLines}]

,

SkipLines=IntegerPart[NoLatticePts/15]+1;

OneString="";

Do[
OneString=StringJoin[OneString,"    1"]
,{j,1,15}];

Do[
WriteString[strm,OneString<>"\r\n"];
,{i,1,SkipLines-1}];

ResidueOnes=Mod[NoLatticePts,15];

OneString="";

Do[
OneString=StringJoin[OneString,"    1"]
,{j,1,ResidueOnes}];

WriteString[strm,OneString<>"\r\n"];

];

If[ValidSOC,

If[YesMessage==True,
Print["Valid SOC matrix provided. Doubling the number of orbitals and writing the Wannier90 data."];
];

RowIndx=0;

Do[
CurrLattPt=CurrRow[[1;;SpaceDim]];

(*Iterating over the data for the given lattice point, interpolating it and loading it into a matrix. We shall double the matrix and, when CurrLattPt is the origin, we shall add to it the SOC matrix.*)
If[ValidSym,
HoppMat=Simplify[ReplaceAll[CurrRow[[SpaceDim+1]],{Symbol[ParamSymbol]->ParamVal}]//N];
,
HoppMat=CurrRow[[SpaceDim+1]];
];

If[CurrLattPt==ConstantArray[0,SpaceDim],
HoppMat2=ZRenorm(ArrayFlatten[IdentityMatrix[2]\[TensorProduct]HoppMat]+SOC)- DiagonalMatrix[ConstantArray[EFermi,2HamltnDim]];
,
HoppMat2=ZRenorm ArrayFlatten[IdentityMatrix[2]\[TensorProduct]HoppMat];
];

Do[
Do[

Do[
WriteString[strm,StringPadLeft[ToString[CurrLattPt[[n]]],5]];
,{n,1,SpaceDim}];

WriteString[strm,StringPadLeft[ToString[l],5],StringPadLeft[ToString[j],5],StringPadLeft[ToString[DecimalForm[Chop[Re[HoppMat2[[l,j]]],Resln]//N,{Length[IntegerDigits[IntegerPart[Re[HoppMat2[[l,j]]]]]]+NoDigits,NoDigits}]],4+1+Length[IntegerDigits[IntegerPart[Re[HoppMat2[[l,j]]]]]]+NoDigits],StringPadLeft[ToString[DecimalForm[Chop[Im[HoppMat2[[l,j]]],Resln]//N,{Length[IntegerDigits[IntegerPart[Im[HoppMat2[[l,j]]]]]]+NoDigits,NoDigits}]],4+1+Length[IntegerDigits[IntegerPart[Re[HoppMat2[[l,j]]]]]]+NoDigits],"\r\n"]

,{l,1,2HamltnDim}]
,{j,1,2HamltnDim}];

,{CurrRow,TBData}];

,

If[YesMessage==True,
Print["Invalid or no SOC matrix provided! Writing an SOC-free Hamiltonian."];
];

Do[
CurrLattPt=CurrRow[[1;;SpaceDim]];

(*Iterating over the data for the given lattice point, interpolating it and loading it into a matrix. Since no valid SOC was provided, we do not double the Hamiltonian matrix.*)
If[ValidSym,
HoppMat=Simplify[ReplaceAll[CurrRow[[SpaceDim+1]],{Symbol[ParamSymbol]->ParamVal}]//N];
,
HoppMat=CurrRow[[SpaceDim+1]];
];

If[CurrLattPt==ConstantArray[0,SpaceDim],
HoppMat2=ZRenorm HoppMat- DiagonalMatrix[ConstantArray[EFermi,HamltnDim]];
,
HoppMat2=ZRenorm  HoppMat;
];

Do[
Do[

Do[
WriteString[strm,StringPadLeft[ToString[CurrLattPt[[n]]],5]];
,{n,1,SpaceDim}];

WriteString[strm,StringPadLeft[ToString[l],5],StringPadLeft[ToString[j],5],StringPadLeft[ToString[DecimalForm[Chop[Re[HoppMat2[[l,j]]],Resln]//N,{Length[IntegerDigits[IntegerPart[Re[HoppMat2[[l,j]]]]]]+NoDigits,NoDigits}]],4+1+Length[IntegerDigits[IntegerPart[Re[HoppMat2[[l,j]]]]]]+NoDigits],StringPadLeft[ToString[DecimalForm[Chop[Im[HoppMat2[[l,j]]],Resln]//N,{Length[IntegerDigits[IntegerPart[Im[HoppMat2[[l,j]]]]]]+NoDigits,NoDigits}]],4+1+Length[IntegerDigits[IntegerPart[Re[HoppMat2[[l,j]]]]]]+NoDigits],"\r\n"]

,{l,1,HamltnDim}]
,{j,1,HamltnDim}];

,{CurrRow,TBData}];
];

If[!YesMessage,
Close[strm];
,
Close[strm]
]
,

If[YesMessage,
Print[Style["Invalid inputs! Try WriteInterpolatedWannier::usage for the specification of correct usage.",Red]]
]
]
]


(* ::Section::Closed:: *)
(*Computing the directional derivative of a band*)


Options[DirectionalDerivative]={"Band"->-1,"Messages"->True,"SanityChecks"->True,"Resolution"->10^-10}

DirectionalDerivative[Hamltn_,kpt_,kdir_,OptionsPattern[]]:=
Module[
{YesMessage,AllBands,\[Lambda]sym,DirDerivs,HDeriv,Degen,DegenFull,Resln,Eigs,EigVals,EigVecs,SpaceDim,HamltnDim},

YesMessage=OptionValue["Messages"];
Resln=10^-10;

If[!BooleanQ[OptionValue["SanityChecks"]] || (BooleanQ[OptionValue["SanityChecks"]] && OptionValue["SanityChecks"]==True),

If[!BooleanQ[OptionValue["Messages"]],
YesMessage=True;
];

If[NumberQ[OptionValue["Resolution"]] && 0 < OptionValue["Resolution"] <1,
Resln=OptionValue["Resolution"];
];

If[!ListQ[kpt] || !AllTrue[kpt,# \[Element] Reals&],
If[YesMessage,
Print[Style["Invalid k-point!",Red]];
];

Return[]
];

If[!ListQ[kdir] || !AllTrue[kdir,# \[Element] Reals&] || Length[kdir]!=Length[kpt],
If[YesMessage,
Print[Style["Invalid direction vector!",Red]];
];

Return[]
];

If[!HermitianMatrixQ[Hamltn[kpt]] || !AllTrue[Flatten[Hamltn[kpt]//N],NumberQ[#]&],
If[YesMessage,
Print[Style["Invalid Hamiltonian!",Red]];
];

Return[]
];

];

SpaceDim=Length[kpt];

HamltnDim=Length[Hamltn[kpt]];

If[IntegerQ[OptionValue["Band"]] && (0<OptionValue["Band"]<=HamltnDim),
AllBands=False;
,
If[YesMessage,
Print[Style["Invalid band number! Computing the directional derivative for all the bands.",Red]];
];
AllBands=True;
];


Eigs=Eigensystem[0.5 (N[Hamltn[kpt]]+Transpose[Conjugate[N[Hamltn[kpt]]]]),Method->"Direct"];

EigVals=Eigs[[1]];
EigVecs=Eigs[[2]];
(* Separating the eigenvalues and eigenvectors in to separate lists *)

Eigs=SortBy[Table[Append[{EigVals[[i]]},EigVecs[[i]]],{i,1,HamltnDim}],First]; (* Sorting the eigensystem by the eigenvalues, rather than by their magnitudes, which is what Mathematica does by default. *)

EigVals=Eigs[[All,1]];
EigVecs=Eigs[[All,2]];

(*We will organise the bands into groups of degenerate bands*)
DegenFull={{1}}; 

Do[
(*Since the bands have already been sorted, the i^th band has to be compared only to the (i-1)^th band for degeneracy*)
If[Abs[EigVals[[i]] - EigVals[[i-1]]]<=Resln,
DegenFull=Insert[DegenFull,i,{-1,-1}];
,
DegenFull=Insert[DegenFull,{i},-1];
];
,{i,2,HamltnDim}];

If[AllBands,
Degen=DegenFull;
,
Degen={DegenFull[[Position[DegenFull,OptionValue["Band"]][[1,1]]]]};
];

(*We can now compute the directional derivative \[PartialD]H along kdir. Its eigenvalues give us the derivatives of the eigenvalues.*)
HDeriv=Conjugate[EigVecs] . Simplify[ReplaceAll[D[Hamltn[kpt+\[Lambda]sym kdir],\[Lambda]sym],{\[Lambda]sym->0}]] . Transpose[EigVecs];
(*As before, we have to use this above form for E^\[Dagger].\[PartialD]H.E since the eigenvector matrix is stored row wise in Mathematica but matrix multiplications are done the normal way whenever possible.*)
(*We are effectively computing \[PartialD]H in the eigenbasis at kpt*)

DirDerivs=First@Last@Reap[Do[
Sow[Sort[Eigenvalues[0.5(HDeriv[[bds,bds]]+ConjugateTranspose[HDeriv[[bds,bds]]])]]];
,{bds,Degen}]];

Chop[DirDerivs,Resln]
]


(* ::Section::Closed:: *)
(*Finding the critical point on a line connecting two given points*)


Options[LocateCriticalPoint]={"MaxIterations"->12,"Messages"->True,"SanityChecks"->True,"Resolution"->10^-4}

LocateCriticalPoint[Hamltn_,kStart_,kEnd_,BandNo_,OptionsPattern[]]:=
Module[
{klow,khigh,kmid,kdir,dElow,dEmid,dEhigh,Degen,CurrItern,MaxIterns,YesMessage,SpaceDim,HamltnDim,Eigs,EigVecs,EigVals,Resln,PtFound},
(*We will use binary search to find a critical point located on the straight line connecting klow and khigh*)

YesMessage=True;
Resln=10^-4;
MaxIterns=12;
klow=kStart;
khigh=kEnd;

(*To save time during the repeated application of this routine, the user may decide to not perform these sanity checks.*)
If[!BooleanQ["SanityChecks"] || (BooleanQ[OptionValue["SanityChecks"]] && OptionValue["SanityChecks"]==True),

If[BooleanQ[OptionValue["Messages"]],
YesMessage=OptionValue["Messages"];
];

If[IntegerQ[OptionValue["MaxIterations"]] && OptionValue["MaxIterations"] > 0,
MaxIterns=OptionValue["MaxIterations"];
If[YesMessage,
Print["The value of MaxIterations has been reset to ",MaxIterns,"."];
]
];

If[NumberQ[OptionValue["Resolution"]] && 0 < OptionValue["Resolution"] <1,
Resln=OptionValue["Resolution"];
,
Resln=10^-10;
];

If[!ListQ[klow] || !AllTrue[klow,# \[Element] Reals&],
If[YesMessage,
Print[Style["Invalid starting point of the range!",Red]];
];

Return[]
];

SpaceDim=Length[klow];

If[!ListQ[khigh] || !AllTrue[khigh,# \[Element] Reals&] || Length[khigh]!=Length[klow],
If[YesMessage,
Print[Style["Invalid ending point of the range!",Red]];
];

Return[]
];

If[!HermitianMatrixQ[Hamltn[ConstantArray[0,SpaceDim]]] || !AllTrue[Flatten[Hamltn[ConstantArray[0,SpaceDim]]//N],NumberQ[#]&],
If[YesMessage,
Print[Style["Invalid Hamiltonian!",Red]];
];

Return[]
];

HamltnDim=Length[Hamltn[ConstantArray[0,SpaceDim]]];

If[!IntegerQ[BandNo] || !(0<BandNo<=HamltnDim),
If[YesMessage,
Print[Style["Invalid band number!",Red]];
];

Return[]
];
];

(*We can now commence the Binary Search*)
(*In order to compute the derivative of the correct band in case of degeneracies, we have to identify the degenerate set first.*)

kdir=Normalize[khigh-klow];
(*We first look at klow*)
dElow=DirectionalDerivative[Hamltn,klow,kdir,"Band"->BandNo,"Messages"->YesMessage,"SanityChecks"->False][[1,1]];
dEhigh=DirectionalDerivative[Hamltn,khigh,kdir,"Band"->BandNo,"Messages"->YesMessage,"SanityChecks"->False][[1,1]];
PtFound=False;
CurrItern=0;
(*
Print[dElow];
Print[dEhigh];*)

If[Chop[dElow,Resln]==0,
Return[{klow,True}]
,
If[Chop[dEhigh,Resln]==0,
Return[{khigh,True}]
,

If[Sign[dEhigh]==Sign[dElow],
If[YesMessage,
Print[Style["The sign of the derivative is the same at the starting and ending points! We can not perform a binary search to locate a critical point in between.",Red]];
];
Return[{{},False}];
,

While[!PtFound && CurrItern<MaxIterns,
kmid=(klow+khigh)/2;
dEmid=DirectionalDerivative[Hamltn,kmid,kdir,"Band"->BandNo,"Messages"->YesMessage,"SanityChecks"->False][[1,1]];
CurrItern=CurrItern+1;

If[Chop[dEmid,Resln]==0,
If[YesMessage,
Print["Point found after ",CurrItern," iteration(s)."]
];
PtFound=True;
,
If[Sign[dEmid]==Sign[dElow],
klow=kmid;
dElow=dEmid;
,
khigh=kmid;
dEhigh=dEmid;
]
]

]

]
]
];

If[PtFound,
Return[{kmid,True}]
,
If[YesMessage,
Print["Could not locate critical point after ",CurrItern," iterations. The value of the derivatives after the last iteration were computed to be ",dElow,", ",dEmid,", ",dEhigh,"."];
];
Return[{{},False}]]
]
