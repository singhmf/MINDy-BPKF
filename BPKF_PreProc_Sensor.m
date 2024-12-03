function[Ymeas,Out,OutStruct]=BPKF_PreProc_Sensor(Ymeas,SourceMesh,ElecPos3d,ProcStruct)
%% Note--this function doesn't rescale or re-center Ymeas

%% ElecPos should be positions on the mesh (or very near)

%% Get nearest point from electrode to mesh
if isfield(SourceMesh,'bnd')
    SourceMesh=SourceMesh.bnd;
end

%% Parcellate & Project data + estimate noise covariances
arQexp=ProcStruct.arQ;
arRexp=ProcStruct.arR;
arHexp=ProcStruct.arH;
%% varQ also gives the additional local contribution to variance (diagonal of vertex-wise Q)
varQ=ProcStruct.varQ;
varR=ProcStruct.varR;

%% arscale scales the full auto-regressive component + additional diagonal (varQ and varR) 
arQscale=ProcStruct.arQscale;
arRscale=ProcStruct.arRscale;


EigThresh=ProcStruct.EigThresh;


if size(ElecPos3d,2)~=3 && size(ElecPos3d,1)==3
    disp('Flipping ElecPos')
    ElecPos3d=ElecPos3d';
end



nChanOrig=size(ElecPos3d,1);

Elec2Vtx=zeros(nChanOrig,1);

%% Find closest source vertex to electrode
for ii=1:nChanOrig
    [~,Elec2Vtx(ii)]=min(sum((SourceMesh.pos-ElecPos3d(ii,:)).^2,2));
end

if numel(unique(Elec2Vtx))~=nChanOrig
    error('Overlapping electrode projections')
end

minDist=BPKF_Mesh_Distances(SourceMesh,ProcStruct.doReScale);

minDist=minDist(Elec2Vtx,Elec2Vtx);


if arHexp~=0
Hset=exp(-arHexp*(minDist.^ProcStruct.orderH));
else
Hset=eye(nChanOrig);
end

if arRexp~=0
Rset=exp(-arRexp*(minDist.^ProcStruct.orderR))+varR*eye(nChanOrig);
else
Rset=varR*eye(nChanOrig);
end
if arQexp~=0
parcQ=arQscale*(exp(-arQexp*(minDist.^ProcStruct.orderQ))+varQ*eye(nChanOrig));
else
parcQ=varQ*eye(nChanOrig);
end


    if ~iscell(Ymeas)
        Ymeas={Ymeas};
    end

if strcmpi(ProcStruct.doMeanRef,'y')
    reRef=eye(nChanOrig)-ones(nChanOrig)/nChanOrig;
    [u,~]=eig(reRef);
    reRef=u(:,2:end)';
    Rset=reRef*Rset*reRef';
    Hset=reRef*Hset;
    for ii=1:numel(Ymeas)
    Ymeas{ii}=reRef*Ymeas{ii};
    end
    Out.reRef=reRef;
end

[Ymeas,OutStruct]=BPKF_PreProc_Proj(Ymeas,Hset,EigThresh,Rset);


if ~iscell(OutStruct.dropCov)
nDimDrop=nChanOrig-size(OutStruct.Hnew,1);

%% Dropped dimensions of measurement are pure noise--estimate original variance
ProjvarR=trace(OutStruct.dropCov)/nDimDrop;

Rmat=OutStruct.Rnew*arRscale+ProjvarR*OutStruct.ProjMeas*OutStruct.ProjMeas';
else
    ProjvarR=cell(1,numel(OutStruct.dropCov));
    Rmat=cell(1,numel(OutStruct.dropCov));
    for iRun=1:numel(OutStruct.dropCov)
nDimDrop=nChanOrig-size(BPKF_Neutral2cell(OutStruct.Hnew,iRun),1);

%% Dropped dimensions of measurement are pure noise--estimate original variance
ProjvarR{iRun}=trace(OutStruct.dropCov{iRun})/nDimDrop;
Rmat{iRun}=BPKF_Neutral2cell(OutStruct.Rnew,iRun)*arRscale+ProjvarR{iRun}*BPKF_Neutral2cell(OutStruct.ProjMeas,iRun)*BPKF_Neutral2cell(OutStruct.ProjMeas,iRun)';
    end
end
Out.H=OutStruct.Hnew;
Out.Q=parcQ;
Out.R=Rmat;
end


