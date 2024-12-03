function[Out]=BPKF_Discrete_Sample(nPt,probVec)
probVec=reshape(probVec,1,numel(probVec));
cumProb=cumsum(probVec);
if cumProb(end)<.999 || cumProb(end)>1.001
    error('Distribution does not add to 1')
end
Out=zeros(1,nPt);
samp=rand(1,nPt);
cumProb=[0 cumProb(1:end-1)];
for ii=find(probVec~=0)
    Out(cumProb(ii)<=samp)=ii;
end
end