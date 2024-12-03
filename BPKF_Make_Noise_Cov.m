function[Qvert,minDist]=BPKF_Make_Noise_Cov(SourceMesh,qExpC,doReScale,arOrder)
%% SourceMesh contains .pos,.tri
%% qExpC (cell or array) contains exponents for kernels
%% doReScale (y/n) denotes whether to rescale distances so max range=1 in largest directions or use native units.
%% arOrder=1: exponential    arOrder=2: gaussian
if nargin==3
    arOrder=1;
end

qSize=size(qExpC);
if iscell(qExpC)
    qExpC=[qExpC{:}]';
    if numel(qExpC)==1
        was1cell=true;
    else
        was1cell=false;
    end
else
    was1cell=false;
end


minDist=BPKF_Mesh_Distances(SourceMesh,doReScale);


if arOrder==2
    minDist=minDist.^2;
end

if numel(qExpC)==1 && ~was1cell
    Qvert=exp(-minDist*qExpC);
else
    Qvert=cell(qSize);
    for iExp=1:numel(qExpC)
Qvert{iExp}=exp(-minDist*qExpC(iExp));
    end
end
end
