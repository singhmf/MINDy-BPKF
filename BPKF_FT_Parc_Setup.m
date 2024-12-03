function[parcLF,parcQ,vecLF,parcellationMat]=BPKF_FT_Parc_Setup(LF,SourceMesh,nYeo,Qvertex)
%% LF=Leadfield (1 x nVert cell)(nChanx3)
%% SourceMesh contains .pos and .tri
%% nYeo=number of Schaeffer parcels
%% Qvertex=process-noise covariance at vertex level

%% parcLF=parcellated leadfield (nChan x nParc)
%% parcQ=parcellated process covariance (nParc x nParc)


%% Directory with Schaeffer parcellation for fs4k, fs5k

parcDir='C:\Users\Matthew\Desktop\HCP\Yeo_Mat';
if size(SourceMesh.pos,1)==20484
    suffRes='_5k.mat';
    disp('Using 5k resolution')
else
    suffRes='_4k.mat';
end

if str2double(nYeo)<100
    nYeo=100*nYeo;
end
combParc=load(fullfile(parcDir,strcat('Yeo_',num2str(nYeo),suffRes)));
combParc=[combParc.yeoL;combParc.yeoR];


%% Calculate source normals
FTnorm=FT_normals(SourceMesh.pos,SourceMesh.tri,'vertex');

%% Multiply each source location by its normal+ rescale for area
vecLF=zeros(size(LF{1},1),size(FTnorm,1));
for ii=1:numel(LF)
    if ~isempty(LF{ii})
    vecLF(:,ii)=LF{ii}*(FTnorm(ii,:)');
    end
end
parcellationMat=(combParc==(1:nYeo))';
parcLF=vecLF*parcellationMat';
if nargin==4
    parcQ=parcellationMat*Qvertex*parcellationMat';
end
end

