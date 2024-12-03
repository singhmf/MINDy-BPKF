function[OutU,OutC,OutX]=BPKF_Separate_Bands(Xdat,Fs,varargin)
bandSet=varargin;
dropPart=1000;
OutU=cell(1,numel(bandSet));
OutC=cell(1,numel(bandSet));
OutX=cell(1,numel(bandSet));
if size(Xdat,2)<(3*dropPart)
    error('not enough data (change drop length')
end
for iB=1:numel(bandSet)
    tmp=bandpass(Xdat',bandSet{iB},Fs);
    tmp=tmp(dropPart:(end-dropPart),:);
    OutX{iB}=tmp';
    tmpC=cov(tmp);
    OutC{iB}=tmpC;
    if nargout==2
        OutC{iB}=tmpC;
    end
    [u,~]=eig(tmpC);
    OutU{iB}=u(:,end-1:end);
end