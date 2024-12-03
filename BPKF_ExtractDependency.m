function[W,diagD,C,SS,V,rtQ,BX,BY]=BPKF_ExtractDependency(ParamValDep)

%if ~isempty(ParamValDep.DepStr)||~any(ParamValDep.ParamStr.fixedVec,[1 2])
ParamStr=BPKF_Apply_Dependency(ParamValDep);


W=ParamStr.W;
diagD=ParamStr.D;
C=ParamStr.C;
SS=ParamStr.S;
V=ParamStr.V;
if ParamStr.doInput
    BX=ParamStr.BX;
    BY=ParamStr.BY;
else
    BX=[];
    BY=[];
end
if ParamStr.isFixed.Qchol
    rtQ=ParamStr.cholQorig;
else
    rtQ=diag(ParamStr.Qchol);
end
end