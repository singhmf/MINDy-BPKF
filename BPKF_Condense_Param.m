function[Param,rtQ]=BPKF_Condense_Param(ParamValDep,doInput)
if doInput
[A,diagD,C,SS,V,rtQ,BX,BY]=BPKF_ExtractDependency(ParamValDep);
else
[A,diagD,C,SS,V,rtQ]=BPKF_ExtractDependency(ParamValDep);
end
if doInput
Param={A,diagD,C,SS,V,BX,BY};
else
    Param={A,diagD,C,SS,V};
end
end