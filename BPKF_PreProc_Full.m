function[Ymeas,Out,OutStruct]=BPKF_PreProc_Full(Ymeas,SourceMesh,LFvertex3d,nYeo,ProcStruct)
%% Reference Matrix
RefMatrix=ProcStruct.RefMatrix;

%% Parcellate & Project data + estimate noise covariances
arQexp=ProcStruct.arQ;
arRexp=ProcStruct.arR;

%% varQ also gives the additional local contribution to variance (diagonal of vertex-wise Q)
varQ=ProcStruct.varQ;

arQscale=ProcStruct.arQscale;
arRscale=ProcStruct.arRscale;

EigThresh=ProcStruct.EigThresh;
%% Whether to rescale based upon "spatial area" per vertex (all connected surfaces)
doSpatialNorm=strcmpi(ProcStruct.doSpatialNorm(1),'y');

if EigThresh>1
    disp('EigThresh>1 assuming you wanted to invert')
    EigThresh=1/EigThresh;
end

LFempty=cellfun(@isempty,LFvertex3d);

%% Get geodesic (surface) convolution of spatial-AR kernels
[expMat,~]=BPKF_Make_Noise_Cov(SourceMesh,[arQexp,arRexp],'y',2);


Qvertex=expMat{1}+varQ*eye(size(expMat{1}));


Rvertex=expMat{2}(~LFempty,~LFempty);

if doSpatialNorm
    disp('Normalizing contributions by surface area')
%% Calculate surface area of phases and summed area of all faces per vertex 
PosArea=BPKF_Surface_Area(SourceMesh);
PosArea=PosArea/mean(PosArea(~LFempty));
NormArea=sqrt(PosArea).*sqrt(PosArea');
%% Only normalize Q+LF here (not R b/c R will be renormalized later by LF)
Qvertex=Qvertex.*NormArea;
for iP=(find(~LFempty))
    LFvertex3d{iP}=LFvertex3d{iP}.*PosArea(iP);
end

end

%% Autoregressive component of R
arRmat=zeros(size(LFvertex3d{find(~LFempty,1)},1));

%% Assuming spherical random noise orientations
for ii=1:3
tmp=Uncellfun(@(xx)(xx(:,ii)),LFvertex3d(~LFempty));
LFdirection=[tmp{:}];
arRmat=arRmat+LFdirection*Rvertex*LFdirection';
end

[parcLF,parcQ]=BPKF_FT_Parc_Setup(LFvertex3d,SourceMesh,nYeo,Qvertex);
nChanOrig=size(parcLF,1);
parcQ=parcQ*arQscale;

%%
if iscell(RefMatrix)
    arRmat0=arRmat;
    arRmat=cell(size(RefMatrix));
    for iC=1:numel(RefMatrix)
arRmat{iC}=RefMatrix{iC}*arRmat0*RefMatrix{iC}';
    end
else
arRmat=RefMatrix*arRmat*RefMatrix';
end
[Ymeas,OutStruct]=BPKF_PreProc_Proj(Ymeas,parcLF,EigThresh,arRmat);


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
%OutStruct.ProjMeas=ProjMeas;
%OutStruct.Rnew=Rnew;
%OutStruct.Yscale=Yscale;
%OutStruct.dropCov=dropCov;
%OutStruct.origCov=origCov;

end
