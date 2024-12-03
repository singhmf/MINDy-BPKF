function[HyperVec,HyperMark,HyperCellMark]=BPKF_PBT_Vectorize(DataStr,fNames)

if ~iscell(fNames)
    fNames={fNames};
end
HyperCellMark=cell(1,numel(fNames));

HyperVec=cell(1,numel(fNames));

for iF=1:numel(fNames)
    fComb=strsplit(fNames{iF},'.');
    tmpVal=getfield(DataStr,fComb{:});
    if iscell(tmpVal)
        tmpVal=Uncellfun(@(xx)(reshape(xx,1,numel(xx))),tmpVal);
        tmpVec=[tmpVal{:}];
        HyperVec{iF}=tmpVec;
        HyperCellMark{iF}=cellfun(@numel,tmpVal);
    else
    HyperVec{iF}=reshape(tmpVal(:),1,numel(tmpVal));
    end
end
HyperMark=cellfun(@numel,HyperVec);
HyperVec=[HyperVec{:}];
end