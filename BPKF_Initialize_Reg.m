function[ParamStr,Reg,nadamStr]=BPKF_Initialize_Reg(ModelSpec,ParamStr,nadamStr)
%% Setup regularization functions


%% For now assuming same wPC for each MINDy block

isYes=@(xx)(isfield(ModelSpec,xx)&&strcmpi(ModelSpec.(xx),'y'));

ParamStr.doReg=isfield(ModelSpec,'Reg')&&~isempty(ModelSpec.Reg);
rtRand=sqrt(ModelSpec.RandScale);

if ~ParamStr.doReg
    Reg=[];
else
    %RegStr=ModelSpec.Reg;
    Ltypes={'L1','L2'};
    for iL0=1:numel(Ltypes)
        iL=Ltypes{iL0};
        Lname=strcat(iL,'names');
    if ~isfield(ModelSpec.Reg,iL)
        Reg.(Lname)={};
    elseif ~isstruct(ModelSpec.Reg.(iL))
        disp(['Assuming ',iL,' regularization refers to W']);
        Reg.(Lname)={'W'};
        Reg.(iL).coeff={ModelSpec.Reg.(iL)};
        Reg.(iL).mean={0};
    else
        Reg.(iL)=ModelSpec.Reg.(iL);
        if ~isfield(Reg.(iL),'mean')
            Reg.(iL).mean=repmat({0},size(Reg.(iL).coeff));
        end
    end
    end
    ParamStr.Reg=Reg;
end
doMINDy=isYes('MINDy');
Reg.doMINDy=doMINDy;
if doMINDy
Mdat=ModelStr.MINDy;
if ~isfield(Mdat,'BlockSz')
    disp('BlockSz not specified: assuming single (full) MINDy block')
    MINDy.PreBlock=[];
    MINDy.PostBlock=[];
    MINDy.PreGroup=1;
    MINDy.PostGroup=1;
    nX=size(ParamStr.W,1);
    wPC=Mdat.wPC(1);
    if ~isfield(Mdat,'Sign')
        disp('Assuming MINDy blocks are unsigned')
        Mdat.Sign=nan;
    end
    mSign=Mdat.Sign;
    if mSign==0
        disp('Interpreting MINDy sign of zero as unsigned (use nan next time)');
        mSign=nan;
    end
    if isnan(mSign)
    ParamStr.Wk1={randn(nX,wPC)*rtRand};
    ParamStr.Wk2={randn(wPC,nX)*rtRand};
    else
        ParamStr.Wk1={abs(randn(nX,wPC))*mSign*rtRand};
        ParamStr.Wk2={abs(randn(wPC,nX))*rtRand};
    end
else
BlockInd=[0 cumsum(Mdat.BlockSz)];
usedBlock=(Mdat.Pre~=0)+(MDat.Post~=0);
if any(usedBlock(:)==1)
    error('Mismatch in MINDy.Pre and MINDy.Post blocks');
end
usedBlock00=find(usedBlock>0)';
[aTmp,bTmp]=find(usedBlock>0);
tmpBlock=zeros(size(usedBlock));
MINDy.BlockPositions=[aTmp,bTmp];
%MINDy.Sign=zeros(1,numel(usedBlock00));


numPre=max(Mdat.Pre,'all');
numPost=max(Mdat.Post,'all');


%if numel(Mdat.wPC)==1
%    Mdat.wPC=repmat(Mdat.wPC,1,max(numPre,numPost));
%en



for ii=1:numel(usedBlock00)
    MINDy.Block{1}(ii,:)=BlockInd(aTmp(ii)+[0 1])+[1 0];
    MINDy.Block{2}(ii,:)=BlockInd(bTmp(ii)+[0 1])+[1 0];
    MINDy.PreGroup(1,ii)=abs(Mdat.Pre(usedBlock00(ii)));
    MINDy.PostGroup(1,ii)=abs(Mdat.Post(usedBlock00(ii)));
    MINDy.Sign{1,ii}=Mdat.Sign(usedBlock00(ii));
end


MINDy.wPC=Mdat.wPC;
if numel(unique(Mdat.Pre(Mdat.Pre~=0)))~=numPre
    error(['Missing MINDy Pre group numbers ',num2str(setdiff(1:numPre,unique(Mdat.Pre(Mdat.Pre~=0))))]);
end
if numel(unique(Mdat.Post(Mdat.Post~=0)))~=numPost
    error(['Missing MINDy Post group numbers ',num2str(setdiff(1:numPost,unique(Mdat.Post(Mdat.Post~=0))))]);
end




for ii=1:numPre
    aTmp=MINDy.PreGroup==ii;
    nPre=unique(MINDy.Block{1}(aTmp,2)-MINDy.Block{1}(aTmp,1));
    %PreInds={MINDy.Block{1}(MINDy.PreGroup==ii,:),MINDy.Block{2}(MINDy.PreGroup==ii)};
    if numel(nPre)~=1
        error(['Inconsistent MINDy Preblock sizes for group ',num2str(ii)]);
    end
    MINDy.nPre(1,ii)=nPre;
    
    if isnan(MINDy.Sign{ii})
    ParamStr.Wk1{ii}=randn(nPre,MINDy.wPC)*rtRand;
    else
        ParamStr.Wk1{ii}=MINDy.Sign{ii}.*abs(randn(nPre,MINDy.wPC)*rtRand);
    end
end
for ii=1:numPost
    aTmp=MINDy.PostGroup==ii;    
    nPost=unique(MINDy.Block{2}(aTmp,2)-MINDy.Block{2}(aTmp,1));
    if numel(nPost)~=1
        error(['Inconsistent MINDy Postblock sizes for group ',num2str(ii)]);
    end
    MINDy.nPost(1,ii)=nPost;
    %% Get signs of all Pre-blocks paired with this post-block
    preSigns=MINDy.Sign{unique(MINDy.PreGroup(aTmp))};
    %% Must be consistent in that either all are nan or all are +/- 1
    if mean(isnan(preSigns),'all')==1
    ParamStr.Wk2{ii}=randn(MINDy.wPC,nPost)*rtRand;
    elseif mean(isnan(preSigns),'all')==0
        ParamStr.Wk2{ii}=abs(randn(MINDy.wPC,nPost))*rtRand;
    else
        error(['MINDy must have either all signed or all unsigned for Preblocks corresponding to Postblock ',num2str(ii)]);
    end
end
end
nadamStr.mWk1=Uncellfun(@(xx)(0*xx),ParamStr.Wk1);
nadamStr.nWk1=nadamStr.mWk1;


nadamStr.mWk2=Uncellfun(@(xx)(0*xx),ParamStr.Wk2);
nadamStr.nWk2=nadamStr.mWk2;
Reg.MINDy=MINDy;
end
end






%function[tmp00]=Tmp_Ensure_Cell(tmp00,sNames)
%if ~iscell(sNames)
%    sNames={sNames};
%end
%for iSname=1:numel(sNames)
%    if ~iscell(tmp00)
%        tmp00.(sNames{iSname})={tmp00.(sNames{iSname})};
%    end
%end
%end