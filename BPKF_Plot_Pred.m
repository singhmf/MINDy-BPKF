function[]=BPKF_Plot_Pred(Out,PredInds)
sz1=floor(sqrt(numel(PredInds)));
sz2=ceil(numel(PredInds)/sz1);
figure
for ii=1:numel(PredInds)
    subplot(sz1,sz2,ii)
    if PredInds(ii)<1
        plot(Out.Pred{end+PredInds(ii)}','r')
        hold on;plot(Out.truePred{end+PredInds(ii)}','.b')
    else
        plot(Out.Pred{PredInds(ii)}','r')
        hold on;plot(Out.truePred{PredInds(ii)}','.b')
    end
end
end