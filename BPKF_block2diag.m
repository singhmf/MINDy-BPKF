function[Out]=BPKF_block2diag(X,nBlock)
%% Assumes blocks are evenly sized
nX=size(X,1);
BlkMat=reshape(1:(nBlock^2),nBlock,nBlock);
OutMark=repelem(BlkMat,nX/nBlock,1);
DiagMark=kron(BlkMat,eye(nX/nBlock));
Out=OutMark;
for ii=1:(nBlock^2)
    Out(OutMark==ii)=X(DiagMark==ii);
end
end
