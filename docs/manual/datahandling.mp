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

Class.G("gInterval")("chrname:string","start:double", "stop:double", "identifier:string")();
Class_noVisibilityMarkers.G;

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

topToBottom(30)(A,B,C,R,P);
S.w = A.e + (50,0);
topToBottom(30)(S,T,I,D,E);
G.w = S.se + (50,0);
topToBottom(20)(G,J,K);

drawObjects(A,B,S,I,T,G,J,K,R,P,D,E,C);
link(inheritance)(B.n -- A.s);
link(aggregation)(S.w -- A.e);
link(aggregation)(I.n -- T.s);
link(aggregation)(pathStepX(S.se, G.nw, 20));
link(aggregation)(rpathStepX(T.e, G.sw, 40));
link(aggregation)(T.n -- S.s);
link(inheritance)(J.n -- G.s);
link(inheritance)(K.n -- J.s);
link(association)(pathStepX(G.sw, R.ne, -20));
link(aggregation)(P.n -- R.s);
link(association)(pathStepX(R.e,D.e, 20));
clink(dashedLink)(D,E);
link(aggregation)(pathStepX(C.e, S.sw, 10));

item(iAssoc)("data")(obj.sw = G.sw + (-20,-10));
item(iAssoc)("belongsTo")(obj.nw = R.e + (20,0));
item(iAssoc)("cdata")(obj.sw = R.ne);
item(iAssoc)("raw")(obj.ne = D.w);

endfig;

bye
