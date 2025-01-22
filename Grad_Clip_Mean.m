function[gradX,startNorm]=Grad_Clip_Mean(gradX,startNorm,iBatch,minBatch,threshScale)
if iscell(gradX)
    totSq=0;
    for ii=1:numel(gradX)
        totSq=totSq+sum(gradX{ii}(:).^2);
    end
else
    totSq=sum(gradX(:).^2);
end


if iBatch<=minBatch
    if iBatch==1
        startNorm=totSq/minBatch;
    else
    startNorm=startNorm+totSq/(minBatch);
    end
else
    Thresh=(startNorm*threshScale);
    if totSq>Thresh
       reScale=sqrt(Thresh/totSq);
    if iscell(gradX)
        for ii=1:numel(gradX)
            gradX{ii}=reScale*gradX{ii};
        end
    else
        gradX=reScale*gradX;
    end
    end
end