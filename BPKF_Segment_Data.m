function[Out,goodOut]=BPKF_Segment_Data(Dat,MarkTimes,winDrop,minLength)
if numel(winDrop)==1
    disp('Assuming symmetric drop')
    winDrop=[-winDrop winDrop];
end
winDrop=reshape(winDrop,1,2);
if winDrop(1)>0||winDrop(2)<0
    disp('Changing drop signs')
    winDrop=[-1 1].*abs(winDrop);
end

if ~iscell(Dat)
    wasCell=false;
    Dat={Dat};
else
    wasCell=true;
end
if ~iscell(MarkTimes)
    MarkTimes={MarkTimes};
end

if numel(MarkTimes)~=numel(Dat)
    error('Incompatible cell dimensions for Dat and Drop')
end


Out=cell(1,numel(Dat));
goodOut=cell(1,numel(Dat));

for iT=1:numel(Dat)
    mTmp=(reshape(MarkTimes{iT},numel(MarkTimes{iT}),1)+winDrop);
    mTmp=unique(max(1,min(mTmp,size(Dat{iT},2))));
    
%    mTmp(mTmp<1)=[];
%    mTmp(mTmp>size(Dat{iT},2))=[];
    mTmp=[0 mTmp' (size(Dat{iT},2)+1)]';
    dTmp=mTmp(2:end)-mTmp(1:(end-1))>minLength;
    goodInd=[1+mTmp(dTmp) mTmp(find(dTmp)+1)-1 ];
    goodOut{iT}=goodInd;
    Out{iT}=cell(1,size(goodInd,1));
    for iK=1:size(goodInd,1)
        Out{iT}{iK}=Dat{iT}(:,goodInd(iK,1):goodInd(iK,2));
    end
end

if ~wasCell
    Out=Out{1};
    goodOut=goodOut{1};
end
end


    

    
