function[ParamStr,nadamStr]=BPKF_Initialize_Param(ModelSpec)

%% REMARK: Can consider at most one BX, BY (not session-specific)

%% TO DO: Add subcortical functionality!

%% Connectivity Coding:
%% +/-1 = full connection
%% +/-2 = diagonal
%% i = unsigned full
%% 2i = unsigned diagonal


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up initial values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RandScale=ModelSpec.RandScale;

usedField=@(xx)(isfield(ModelSpec,xx)&&~isempty(ModelSpec.(xx)));

%% Assuming Q is the last variable name
origVar={'W','D','C','S','V'};
optVar={'Qchol'};

if ModelSpec.doInput
inputVar={'BX','BY'};
nU=ModelSpec.nU;
if iscell(ModelSpec.H)
nH=size(ModelSpec.H{1},1);
else
    nH=size(ModelSpec.H,1);
end
else
    inputVar={};
end
reqVar=[origVar,inputVar];
varNames=[reqVar,optVar];

ParamStr.dropGrad=[];

ParamStr.FixedVal=cell(1,numel(reqVar));
for ii=1:(numel(reqVar))
ParamStr.isFixed.(varNames{ii})=usedField(strcat('fix',varNames{ii}));
if ParamStr.isFixed.(varNames{ii})
    ParamStr.FixedVal{ii}=ModelSpec.(strcat('fix',varNames{ii}));
    if any(isnan(ParamStr.FixedVal{ii}))
        disp(['Interpreting NANs in fix',varNames{ii},' as indicating free-values']);
       ParamStr.FreeMask{ii}=isnan(ParamStr.FixedVal{ii});
       ParamStr.FixedVal{ii}(ParamStr.FreeMask{ii})=0;
    else
        ParamStr.FreeMask{ii}=0;
    end
end
end
ParamStr.isFixed.Qchol=~(isfield(ModelSpec,'freeQ')&&~strcmpi(ModelSpec.freeQ(1),'n'));
if ParamStr.isFixed.Qchol
    ParamStr.dropGrad.Q=0;
end
ParamStr.Qorig=ModelSpec.Q;



%ParamStr.Fixed.W=usedfield('fixW');
%ParamStr.Fixed.C=usedfield('fixC');
%ParamStr.Fixed.D=usedField('fixD');
%ParamStr.Fixed.S=usedField('fixS');
%ParamStr.Fixed.V=usedField('fixV');



if usedField('fixW')
    ParamStr.W=ModelSpec.fixW;
    ParamStr.dropGrad.W=0;
else
if isfield(ModelSpec,'CortConnBlock')&&~isempty(ModelSpec.CortConnBlock)
    ConnBlock=ModelSpec.CortConnBlock;
    Wmask=double(ModelSpec.Wmask);%~=0);
    disp('Generating connectivity via block specification')
    if mean(ismember(double(Wmask(:)),[0 1]))~=1
        disp('Treating Wmask as binary (ignoring sign) as these are specified in the Block Conn');
    end
if any(((imag(ConnBlock)~=0)+(real(ConnBlock)~=0))>1)
    error('ConnBlock entries are ill-specified (nonzero real and imag parts')
end

abR=@(xx)(abs(real(xx)));

nxI=eye(size(Wmask));

FullConn=sign((ConnBlock)).*(abR(ConnBlock)==1);
DiagConn=sign((ConnBlock)).*(abR(ConnBlock)==2);


Conn=kron(FullConn,Wmask)+kron(DiagConn,nxI);
unSgnConn=kron(ConnBlock==1i,Wmask)+kron(ConnBlock==(2i),nxI);

Conn(unSgnConn~=0)=nan;

%% Warning if mixing signs between ConnBlock and W
if any(~ismember(Wmask,[0,1]),[1 2])&&(numel(ConnBlock)>1)
    disp('WARNING: using non-binary Wmask with block-specified connectivity')
end

else
    disp('No block-structure specified, using entered Wmask');
    Conn=ModelSpec.Wmask;
end


ParamStr.TrueMask=Conn;
ParamStr.dropGrad.W=Conn~=0;

sd0=RandScale*randn(size(Conn));
W=NanToZero(Conn).*abs(sd0);
W(isnan(Conn))=sd0(isnan(Conn));
end
ParamStr.W=W;

nX=size(W,1);

ParamStr.C=RandScale*randn(nX,1);
ParamStr.D=min(ModelSpec.maxD,ModelSpec.minD+RandScale*abs(randn(nX,1)));
ParamStr.S=min(ModelSpec.maxS,ModelSpec.minS+RandScale*abs(randn(nX,1)));
ParamStr.V=RandScale*randn(nX,1);


%% Remark: dropGrad is reverse-coded
for ii=1:numel(origVar)
    tmpName=origVar{ii};
    if ~strcmpi(tmpName,'W')
    if ParamStr.isFixed.(tmpName)
        ParamStr.(tmpName)=(ParamStr.FreeMask{ii}.*ParamStr.(tmpName))+ParamStr.FixedVal{ii};
        ParamStr.dropGrad.(tmpName)=ParamStr.FreeMask{ii};
    end
    end
end




if ModelSpec.doInput
    if usedField('fixBX')
    ParamStr.BX=ModelSpec.fixBX;
    ParamStr.dropGrad.BX=ParamStr.FreeMask{strcmpi(ParamStr.varOrder,'BX')};
    else
    bx0=RandScale*randn(nX,nU);
    if usedField('BXmask')
        Bconn=ModelSpec.BXmask;
        ParamStr.dropGrad.BX=Bconn~=0;
BX=NanToZero(Bconn).*abs(bx0);
BX(isnan(Bconn))=bx0(isnan(Bconn));
    else
        BX=bx0;
    end
    ParamStr.BX=BX;
    end
    if usedField('fixBY')
    ParamStr.BY=ModelSpec.fixBY;
    ParamStr.dropGrad.BY=ParamStr.FreeMask{strcmpi(ParamStr.varOrder,'BY')};
    else
    by0=RandScale*randn(nH,nU);
    if usedField('BYmask')
    Bconn=ModelSpec.BYmask;
    ParamStr.dropGrad.BY=Bconn~=0;
    BY=NanToZero(Bconn).*abs(by0);
    BY(isnan(Bconn))=by0(isnan(Bconn));
    else
        BY=by0;
    end
    ParamStr.BY=BY;
    end
end


ParamStr.Qorig=ModelSpec.Q;
ParamStr.cholQorig=chol(ModelSpec.Q,'lower');
if (isfield(ModelSpec,'freeQ')&&~strcmpi(ModelSpec.freeQ(1),'n'))
    ParamStr.Qchol=min(ModelSpec.maxQchol,randn(nX,1)*RandScale+ModelSpec.minQchol);
end
    %else
%    ParamStr.Q=ModelSpec.Q;
%end



fixedVec=false(size(varNames));
for ii=1:numel(varNames)
    fixedVec(ii)=ParamStr.isFixed.(varNames{ii});
end

varText=strcat({' '},varNames,{' '});
%fixedVec=[ParamStr.isFixed.W ParamStr.isFixed.D ParamStr.isFixed.C, ParamStr.isFixed.S, ParamStr.isFixed.V,ParamStr.Fixed.Q];
disp(['Using fixed ',[varText{fixedVec}]]);
ParamStr.fixedVec=fixedVec;
if ~ParamStr.isFixed.Qchol%(isfield(ModelSpec,'FreeQ')&&~strcmpi(ModelSpec.FreeQ(1),'n'))
ParamStr.varOrder=varNames;
else
ParamStr.varOrder=varNames(~strncmpi(varNames,'Q',1));%{'W','D','C','S','V'};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save parameter sign assignments (C,V are unsigned)

limWmax=inf(nX);
limWmin=-inf(nX);
limWmin(Conn==1)=0;
limWmax(Conn==-1)=0;

limWmin(Conn==0)=0;
limWmax(Conn==0)=0;

ParamStr.doMax=false(size(ParamStr.varOrder));
ParamStr.doMin=false(size(ParamStr.varOrder));
for ii=1:numel(ParamStr.varOrder)
    if strcmpi(ParamStr.varOrder{ii},'W')
        ParamStr.doMin(ii)=true;
        ParamStr.doMax(ii)=true;
        ParamStr.minBound{ii}=limWmin;
        ParamStr.maxBound{ii}=limWmax;
    else
if isfield(ModelSpec,strcat('max',ParamStr.varOrder{ii}))
    if ~isempty(ModelSpec.(strcat('max',ParamStr.varOrder{ii})))
    ParamStr.doMax(ii)=true;
    ParamStr.maxBound{ii}=ModelSpec.(strcat('max',ParamStr.varOrder{ii}));
    end
end

if isfield(ModelSpec,strcat('min',ParamStr.varOrder{ii}))
    if ~isempty(ModelSpec.(strcat('min',ParamStr.varOrder{ii})))
    ParamStr.doMin(ii)=true;
    ParamStr.minBound{ii}=ModelSpec.(strcat('min',ParamStr.varOrder{ii}));
    end
end
    end
end


    
%minBound=zeros(1,numel(ParamStr.varOrder)-1);
%maxBound=zeros(1,numel(ParamStr.varOrder)-1);


%for ii=2:numel(minBound)
%    if isfield(ModelSpec,strcat('min',ParamStr.varOrder{ii}))
%        minBound(ii-1)=ModelSpec.(strcat('min',ParamStr.varOrder{ii}));
%    else
%        minBound(ii-1)=-inf;
%    end
%    if isfield(ModelSpec,strcat('max',ParamStr.varOrder{ii}))

%        maxBound(ii-1)=ModelSpec.(strcat('max',ParamStr.varOrder{ii}));
%    else
%        maxBound(ii-1)=inf;
%    end
%end

%signTmp=sign(maxBound);signTmp(signTmp~=sign(minBound))=nan;

%ParamStr.ParamSign=[ParamStr.TrueMask repelem(signTmp,nX,1)];
%ParamStr.minVal=[limWmin repmat(minBound,nX,1)];
%ParamStr.maxVal=[limWmax repmat(maxBound,nX,1)];


%if ~ParamStr.isFixed.Qchol
%    ParamStr.ParamSign=[ParamStr.TrueMask repelem([1 nan 1 nan 1],nX,1)];
%    ParamStr.minVal=[limWmin repmat([ModelSpec.minD -inf ModelSpec.minS -inf ModelSpec.Qchol],nX,1)];
%    ParamStr.maxVal=[limWmax repmat([ModelSpec.maxD inf ModelSpec.maxS inf ModelSpec.Qchol],nX,1)];
%else
%    ParamStr.ParamSign=[ParamStr.TrueMask repelem([1 nan 1 nan],nX,1)];
%    ParamStr.minVal=[limWmin repmat([ ModelSpec.minD -inf ModelSpec.minS -inf],nX,1)];
%    ParamStr.maxVal=[limWmax repmat([ ModelSpec.maxD inf ModelSpec.maxS inf],nX,1)];
%end



for ii=1:numel(ParamStr.varOrder)
    nadamStr.(strcat('m',ParamStr.varOrder{ii}))=zeros(size(ParamStr.(ParamStr.varOrder{ii})));
    nadamStr.(strcat('n',ParamStr.varOrder{ii}))=zeros(size(ParamStr.(ParamStr.varOrder{ii})));
end




end