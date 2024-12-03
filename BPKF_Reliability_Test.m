function[Out]=BPKF_Reliability_Test(ooP1,ooP2,BlockSize)

%% BlockSize is a vector (e.g. [100 100]) denoting size of blocks to do separate
%% reliability over


nX=size(ooP1.Param{1},1);
if nargin==2
    BlockSize=nX;
end
if numel(BlockSize)==1
    BlockSize=repmat(BlockSize,1,nX/BlockSize);
end
BlockInds=cell(size(BlockSize));
BlockInds{1}=1:BlockSize(1);
nCumBlock=cumsum(BlockSize);
for iBlock=2:numel(BlockSize)
    BlockInds{iBlock}=(1:BlockSize(iBlock))+nCumBlock(iBlock-1);
end
nBlock=numel(BlockSize);
Out.BlockInds=BlockInds;
Out.Notes='For matrix parameters, reliability is {offdiag, diag, full}';

Out.isMatrix=false(size(ooP1.Param));
for iP=1:5%numel(ooP1.Param)
    if max(size(ooP1.Param{iP}))>1
    if min(size(ooP1.Param{iP}))~=1
        Out.isMatrix(iP)=true;
        for iBlock=1:nBlock
            for jBlock=1:nBlock
                tmp1=ooP1.Param{iP}(BlockInds{iBlock},BlockInds{jBlock});
                tmp2=ooP2.Param{iP}(BlockInds{iBlock},BlockInds{jBlock});
                nonZero=and(tmp1~=0,tmp2~=0);
                tyMark={NoDiag(nonZero)==1,eye(size(tmp1,1))==1,nonZero};
                for iTy=1:3
                [Out.r{iP}{iTy}(iBlock,jBlock),Out.rho{iP}{iTy}(iBlock,jBlock)]=...
                    dropCorr(tmp1(tyMark{iTy}),tmp2(tyMark{iTy}));
                end
            end
        end
    else
        for iBlock=1:nBlock
                tmp1=ooP1.Param{iP}(BlockInds{iBlock});
                tmp2=ooP2.Param{iP}(BlockInds{iBlock});
                                [Out.r{iP}(iBlock),Out.rho{iP}(iBlock)]=dropCorr(tmp1,tmp2);
        end
    end
    end
end
end




function[rOut,rhoOut]=dropCorr(dat1,dat2)
if isempty(dat1)||isempty(dat2)
    rOut=nan;rhoOut=nan;
else
    rOut=corr(dat1(:),dat2(:));
    rhoOut=corr(dat1(:),dat2(:),'type','spearman');
end
end