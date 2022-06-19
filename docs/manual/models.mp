input metauml;
outputtemplate:="%j%c.mps";
iClass.iMethodStack.top :=0;
iClass.iMethodStack.bottom :=0;

beginfig(1);
Class.P("perStrandInfo")("forward:double", "reverse:double")();
Class_noVisibilityMarkers.P;

Class.R("Responsibilities")("Ri:perStrandInfo", "Rk:perStrandInfo")();
Class_noVisibilityMarkers.R;

Class.N("Normal")("mu:double","sigma:double")();
Class_noVisibilityMarkers.N;

Class.E("Exponential")("lambda:double")();
Class_noVisibilityMarkers.E;

Class.U("Uniform")("lower:double","upper:double")();
Class_noVisibilityMarkers.U;

Class.B("BasicModel")("pi:double", "weight:double", "sufficiency_stats:Responsibilities", "betaPi:Priors", "DirichletW:Priors")();
Class_noVisibilityMarkers.B;

Class.F("FullModel")("bidir:Bidirectional", "forwardElongation:UniformModel", "reverseElongation:UniformModel")();
Class_noVisibilityMarkers.F;

Class.M("Bidirectional")("loading:Normal", "initiation:Exponential", "footprint:double", "footprint:double", "GammaSigma:Priors", "GammaLambda:Priors", "sumExpected:sumOverN")();
Class_noVisibilityMarkers.M;

Class.L("UniformModel")("uni:Uniform")();
Class_noVisibilityMarkers.L;

Class.O("Priors")("alpha:double", "beta:double")();
Class_noVisibilityMarkers.O;

Class.S("sumExpected")("sumRExpY:double", "sumRExpX:double", "sumRExpX2:double", "sumRXExpY:double")();
Class_noVisibilityMarkers.S;

Class.C("bidirConstraints")("lambdaMin", "lambdaMax", "fp_Min", "fp_Max")();
Class_noVisibilityMarkers.C;

topToBottom(20)(N,E,P,R,U);
M.w = N.e + (50,0);
topToBottom(20)(M,B,L,F);
S.w = M.e + (60,0);
topToBottom(20)(S,C,O);

drawObjects(P,R,N,E,U,B,M,F,L,O,S,C);

link(aggregation)(R.n -- P.s);
link(aggregation)(pathStepX(R.e, B.w, 10));
link(aggregation)(pathStepX(B.e, O.w, 20));

link(inheritance)(M.s -- B.n);
link(inheritance)(L.n -- B.s);
link(aggregation)(rpathStepX(N.e, M.w, 10));
link(aggregation)(rpathStepX(E.e, M.w, 20));

link(aggregation)(pathStepX(U.e, L.w, 20));
link(aggregation)(M.e -- S.w);
link(aggregation)(pathStepX(M.e, C.w, 30));

path bigone;
bigone := F.s -- F.s + (0,-10) -- F.s + (-180,-10) -- M.n + (-180,20) -- M.n + (0,20) -- M.n;
link(aggregation)(bigone);
link(aggregation)(pathStepX(B.e, F.e, 20));


endfig;

bye
