function[Out]=BPKF_Extract_Block(W,BlockInd)
Out=cell(1,size(BlockInd,1));
for ii=1:size(BlockInd,1)
    Out{1,ii}=W(BlockInd(ii,1):BlockInd(ii,2),BlockInd(ii,3):BlockInd(ii,4));
end
end
