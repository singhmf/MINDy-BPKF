function[DataStr,OutVals]=BPKF_PBT_Explore(DataStr,fNames,pertVec)

if ~iscell(fNames)
    fNames={fNames};
end

if nargout==2
    OutVals=cell(numel(fNames),1);
end
for iF=1:numel(fNames)
    fComb=strsplit(fNames{iF},'.');
    tmpVal=getfield(DataStr,fComb{:});
    if iscell(tmpVal)
        tmpVal=Uncellfun(@(xx)(xx.*pertVec(randi(numel(pertVec),size(xx)))),tmpVal);
    else
    tmpVal=tmpVal.*pertVec(randi(numel(pertVec),size(tmpVal)));
    end
    if nargout==2
        OutVals{iF}=tmpVal;
    end
    DataStr=setfield(DataStr,fComb{:},tmpVal);
end
end