function[FullMat,DepGroup]=BPKF_Standard_SubCort_Model_00(Wmask0,HemiMark)
%% EI model with sub-cortex; can't use cortical-block format

tmp=sort(unique(HemiMark));
if tmp(2)==2 && tmp(1)==1
    HemiMark=1-HemiMark;
elseif tmp(1)==-1 && tmp(2)==1
    HemiMark=(1+HemiMark)/2;
end

tmp=sort(unique(HemiMark));
if tmp(1)~=0 || tmp(2)~=1
    error('Entries should be {0,1} or {-1,1} or {1,2}')
end
nX=size(Wmask0,1);
if numel(HemiMark)~=nX
    error('Wmask and Hemi incompatible')
end

%% Using 3 subcort regions/hemisphere with 2 nodes each
nThal=3;

HemiMark=reshape(HemiMark,1,numel(HemiMark));

ThalLocal=[ones(nThal) eye(nThal);eye(nThal) ones(nThal)];

Cort2Thal=[repmat(1-HemiMark,nThal,1);repmat(HemiMark,nThal,1)];

Wfull=[Wmask0 Cort2Thal';Cort2Thal ThalLocal];

%% 2 Hems
FullMat=repmat([Wfull -eye(size(Wfull))],2,1);

%% Dependency Groups
subDep=repmat((1:2*nThal)',1,nX);

%DepGroup=kron([1 2;3 4],eye(nX));
%DepGroup=[DepGroup (4+subDep');(4+2*nThal+subDep) zeros(2*nThal)];

%% Exc  | Exc[1] Thal+[1+(1:2n)] Inh Thal-
%% Thal | Exc[(2n+1)+(1:2n)] Thal+ Inh Thal-
nT=2*nThal;
DepGroup=[  eye(nX) repmat(1+(1:nT),nX,1) (2+nT)*eye(nX) zeros(nX,nT);...
            repmat((2+nT)+(1:nT)',1,nX) zeros(nT,nX+2*nT)];
DepGroup=[DepGroup;(DepGroup+max(DepGroup,[],'all')).*(DepGroup~=0)];
end


%FullMat=[Wcort [Cort2Thal' zeros(nX,2*nThal)];repmat(Cort2Thal,2,1) zeros(2*nThal,nX) ThalBlock];

