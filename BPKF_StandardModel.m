function[ModelSpec]=BPKF_StandardModel(ModelClass,DependencyClass,ModelSpec)

if nargin<2
    ModelSpec=[];
    if nargin==1
        DependencyClass=[];
    end
end

%% Constant dependencies=same value for all elements within group:
%%   block array of indices by subpop with sign or complex (for free sign)
%%   0=no constant dependency

%% Dep Order [W D C S V diagQ]

%% Connectivity Coding:
%% +/-1 = full connection
%% +/-2 = diagonal
%% i = unsigned full
%% 2i = unsigned diagonal

switch ModelClass
    case 0
        disp('Doing Unsigned Model')
        ModelSpec.Hkron=1;
        ModelSpec.CortConnBlock=1i;
    case 1
        disp('Doing E-I Model')
        ModelSpec.Hkron=[1 0];
        ModelSpec.CortConnBlock=[1 -2;1 -2];
    case 2
        disp('Doing E-(E-I) Model')
        ModelSpec.Hkron=[1 1 0];
        ModelSpec.CortConnBlock=[1 2 -2;1 2 -2;2 2 -2];
end
nPop=size(ModelSpec.Hkron,2);
if ~isempty(DependencyClass)&&DependencyClass~=0
    switch DependencyClass
        case 1
            disp('Doing region-specific for recurrents only: All reccurrents ind')
            Wpart=reshape(1:(nPop^2),nPop,nPop);
            tmp=sign(ModelSpec.CortConnBlock);tmp(isnan(tmp))=1j;
            ModelSpec.Group=[Wpart.*tmp zeros(nPop,6)];
            ModelSpec.Const=[zeros(size(Wpart)),reshape(1:(5*nPop),nPop,5)];
        case 2
            disp('Common Subspace')
            ModelSpec.dfGroup=3;
            ModelSpec.Const=[0 0 1 0 2 3];
            ModelSpec.Group=[1 1 0 1 0 0];
    end
end

end