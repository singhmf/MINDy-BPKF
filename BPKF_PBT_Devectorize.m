function[DataStr]=BPKF_PBT_Devectorize(HyperVec,DataStr,fNames)
   

   

if ~iscell(fNames)
    fNames={fNames};
end




for iF=1:numel(fNames)
    fComb=strsplit(fNames{iF},'.');
    tmpVal=getfield(DataStr,fComb{:});
    if iscell(tmpVal)
        for jj=1:numel(tmpVal)
            tmpVal{jj}(:)=HyperVec(1:numel(tmpVal{jj}));
            HyperVec(1:numel(tmpVal{jj}))=[];
        end
    else
        tmpVal(:)=HyperVec(1:numel(tmpVal));
        HyperVec(1:numel(tmpVal))=[];
    end
    DataStr=setfield(DataStr,fComb{:},tmpVal);
end
end