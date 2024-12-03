function[]=BPKF_Plot_Chan_Pred(Out,PredInds)
for ii=1:numel(PredInds)
    figure
    
    if PredInds(ii)<1
        iInd=numel(Out.Pred)+PredInds(ii);
    else
        iInd=PredInds(ii);
    end
        nX=size(Out.Pred{iInd},1);
        
    sz1=floor(sqrt(nX));
    sz2=ceil(nX/sz1);

    for jj=1:nX
        subplot(sz1,sz2,jj)
        plot(Out.Pred{iInd}(jj,:),'r')
        hold on;plot(Out.truePred{iInd}(jj,:),'.b')
    end
end