input metauml;
outputtemplate:="%j%c.mps";
iClass.iMethodStack.top :=0;
iClass.iMethodStack.bottom :=0;

beginfig(1);
Class.W("ModelWrapper")("type:ModTypes", "*gene:FullModel", "*bidir:Bidirectional")();
Class_noVisibilityMarkers.W;

Class.X("ModelContainer")("noise:UniformModel", "K:int", "setModels:vector{ModelWrapper}", "ll:double", "limits:bidirConstraints")();
Class_noVisibilityMarkers.X;

Class.A("AlgorithmControl")("convergence_threshold:double", "max_iterations:int", "seed:bool", "noise_max:double", "r_mu:double", "elon_move:bool", "maxUniformIter:int")();
Class_noVisibilityMarkers.A;

Class.E("EMalg")("converged:bool", "control:AlgorithmControl", "models:ModelContainer")();
Class_noVisibilityMarkers.E;

Note.D("dInterval");

Class.S("EMseeds")("mu_seeds:vector{double}", "weights:vector{double}", "numgen:Random")();
Class_noVisibilityMarkers.S;

topToBottom(20)(W,X);
A.w = W.e + (70,0);
topToBottom(20)(A,E);
S.w = A.e + (30,0);
topToBottom(20)(S,D);

drawObjects(W,X,A,E,D,S);

link(aggregation)(X.n -- W.s);
link(aggregation)(E.n -- A.s);
link(aggregation)(rpathStepX(X.e, E.w, 30));
link(dependency)(D.w -- E.e);
link(association)(S.s -- D.n);

endfig;

bye
