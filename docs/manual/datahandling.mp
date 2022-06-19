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

Class.S("SetROI")("bimap:object","gInterval:map{vector}")();
Class_noVisibilityMarkers.S;

Class.I("ITree")("left:ptr","right:ptr", "midpoint:double", "node:ptr")();
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

topToBottom(20)(A,B,C);
S.w = A.e + (50,0);
topToBottom(20)(S,I);
G.w = S.se + (50,0);
topToBottom(20)(G,J,K);
R.w = G.e + (50,0);
topToBottom(20)(R,P);
D.w = R.e + (50,0);
E.s = D.nw + (10,10);

drawObjects(A,B,S,I,G,J,K,R,P,D,E,C);
link(inheritance)(B.n -- A.s);
link(aggregation)(S.w -- A.e);
link(aggregation)(pathStepX(S.se, G.nw, 20));
link(aggregation)(rpathManhattanX(I.e, G.w));
link(aggregation)(I.n -- S.s);
link(inheritance)(J.n -- G.s);
link(inheritance)(K.n -- J.s);
clink(association)(G, R);
link(aggregation)(P.n -- R.s);
clink(association)(R,D);
clink(dashedLink)(D,E);
link(aggregation)(pathStepX(C.e, S.sw, 10));

item(iAssoc)("data")(obj.sw = G.e);
item(iAssoc)("belongsTo")(obj.ne = R.w);
item(iAssoc)("cdata")(obj.sw = R.e);
item(iAssoc)("raw")(obj.ne = D.w);

endfig;

bye
