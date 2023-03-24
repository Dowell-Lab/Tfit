input metauml;
outputtemplate:="%j%c.mps";
iClass.iMethodStack.top :=0;
iClass.iMethodStack.bottom :=0;

beginfig(1);
Class.A("bedfile")("filename:string","setROI:object")();
Class_noVisibilityMarkers.A;

Class.B("bedGraph")("hasROI:bool")();
Class_noVisibilityMarkers.B;

Class.C("bimap")("num_elements:int", "index2str:map", "str2index:map")();
Class_noVisibilityMarkers.C;

Class.S("SetROI")("bimap:object","gInterval:map<vector{gInterval}>","searchable:map<CITree>", "treesExist:bool")();
Class_noVisibilityMarkers.S;

Class.T("CITree")("root:Inode")();
Class_noVisibilityMarkers.T;

Class.I("ITree")("left:ptr","right:ptr", "midpoint:double", "OverlapCenter:vector{gInterval}")();
Class_noVisibilityMarkers.I;

Class.F("bed4")("chrname:string","start:double", "stop:double", "identifier:string")();
Class_noVisibilityMarkers.F;

Class.J("bed6")("score:double","strand:char")();
Class_noVisibilityMarkers.J;

Class.K("bed12")("Seeds:object")();
Class_noVisibilityMarkers.K;

Class.R("RawData")("PointCov:vector")();
Class_noVisibilityMarkers.R;

Class.P("PointCov")("position:double","coverage:double")();
Class_noVisibilityMarkers.P;

Class.D("dInterval")("X:vector{double}", "scale:double", "delta:double", "bins:int")();
Class_noVisibilityMarkers.D;

Note.E("X : [position, forward_cov, reverse_cov]*bins");

Class.G("gInterval")("raw:RawData", "cdata:dInterval", "gcoords:bed12")();
Class_noVisibilityMarkers.G;

Class.H("EMSeeds")("...")();
Class_noVisibilityMarkers.H;

topToBottom(30)(A,B,C,R,P);
S.w = A.e + (50,0);
topToBottom(30)(S,T,I,D,E);
F.w = S.se + (50,0);
topToBottom(20)(F,J,K,H,G);

drawObjects(A,B,S,I,T,F,J,K,R,P,D,E,C,G,H);
link(inheritance)(B.n -- A.s);
link(aggregation)(H.n -- K.s);
link(aggregation)(S.w -- A.e);
link(aggregation)(I.n -- T.s);
link(aggregation)(pathStepX(S.se, K.nw, 20));
link(aggregationUni)(pathStepX(G.nw, K.sw, -20));
link(aggregation)(pathStepX(G.w, R.ne, -30));
link(aggregation)(pathStepX(G.sw, D.e, -50));
link(aggregation)(T.n -- S.s);
link(inheritance)(J.n -- F.s);
link(inheritance)(K.n -- J.s);
link(aggregation)(P.n -- R.s);
clink(dashedLink)(D,E);
link(aggregation)(pathStepX(C.e, S.sw, 10));

%item(iAssoc)("raw")(obj.ne = D.w);

endfig;

bye
