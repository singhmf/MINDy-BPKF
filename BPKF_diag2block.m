function[Out]=BPKF_diag2block(X,nBlock)
%% input should have width =nBlock
if size(X,2)~=nBlock
    error('Size 2 of input should equal nBlock')
end
nX=size(X,1);
BlkMat=reshape(1:(nBlock^2),nBlock,nBlock);
OutMark=repelem(BlkMat,nX/nBlock,1);
DiagMark=kron(BlkMat,eye(nX/nBlock));
Out=DiagMark;
for ii=1:(nBlock^2)
    Out(DiagMark==ii)=X(OutMark==ii);
end
end
