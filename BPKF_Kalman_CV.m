function[Out]=BPKF_Kalman_CV(ParamValDep,KalSpec,NBatch,MeasSet,Input)

%% cvStr has fields H,Q,R, nStep, nRec

%% 1.4 More variable setup
opts.POSDEF = true;
opts.SYM = true;

Etotal=0;
Pbase=KalSpec.Pbase;
%% 1.7 Kalman Filter Specification
isEKF=strcmpi(KalSpec.KFtype(1),'E');
%isCKF=sum(strcmpi(KalSpec.KFtype(1),{'U','C'}),'all')~=0;
isCKF=sum(strcmpi(KalSpec.KFtype(1),'C'));
isUKF=sum(strcmpi(KalSpec.KFtype(1),'U'));
if ~ismember(KalSpec.KFtype(1),'EC')
    error('KFtype should be EKF, UKF or CKF')
end

if isCKF 
%% Cubature KF
BaseMat=sqrt(nX)*[eye(nX) -eye(nX)];
BaseWeight=repmat(1/(2*nX),1,2*nX);
opts2.LT = true;
InX=eye(nX);

kSgn=[1 -1]/sqrt(nX);
end

if strcmpi(KalSpec.BFcov(1),'s')
    SimLength=KalSpec.SimLength;
    nSim=KalSpec.nSim;
    nSaveStart=KalSpec.nSaveStart;
    simStart=randn(nX,nSaveStart);
    Pfix=KalSpec.Pfix;
    decP=KalSpec.decP;
    decFix=KalSpec.decFix;
    simCovOld=eye(nX);
elseif strcmpi(KalSpec.BFcov(1),'y')
    Pfix=KalSpec.Pfix;
end
    if strcmpi(KalSpec.BFcov(1),'s')
simRand=cell(1,SimLength);
sCovTmp=zeros(nX);
sMeanTmp=zeros(nX,nSim);
    end


if isfield(ModelSpec,'nRep')
    nRep=ModelSpec.nRep;
else
    nRep=1;
end

minStep=KalSpec.minStep;
nStep=KalSpec.nStep;
nRec=KalSpec.nRec;
if isfield(ModelSpec,'dropEnd')
    dropEnd=ModelSpec.dropEnd;
else
    dropEnd=nStep+1+nRec;
end

if nargin==6
    doInput=false;
else
    doInput=isempty(Input);
end

if doInput
[W,diagD,C,SS,V,rtQ,BX,BY]=BPKF_ExtractDependency(ParamValDep);
else
[W,diagD,C,SS,V,rtQ]=BPKF_ExtractDependency(ParamValDep);
end


    SimLength=KalSpec.SimLength;
    nSim=KalSpec.nSim;

    ATS=SS.*(W');
    nX2=2*nX;
    n2V=-2*V;    
    n2S=-2*SS;
    %% 3 Initialize distributions
    if ~strcmpi(KalSpec.BFcov(1),'n')
        %% 3.1 Use simulations (if specified)
    if strcmpi(KalSpec.BFcov(1),'s')
       
    simXC=cell(1,SimLength+1);
    simX=simStart(:,randperm(nSaveStart,nSim));
%    simRand=randn(nX,nSim,SimLength);
    simFunC=cell(1,SimLength);
    sMeanTmp=0*sMeanTmp;
    sCovTmp=0*sCovTmp;
    for iS=1:SimLength
    simRand{iS}=randn(nX,nSim);
    simXC{iS}=simX;
%    SimFun=-1+2./(1+exp(-2*(SS.*simX+V)));
    SimFun=-1+2./(1+exp(n2S.*simX+n2V));
    simFunC{iS}=SimFun;
    if isDiagQ
        simX=W*SimFun+diagD.*simX+C+dgRtQ.*simRand{iS};%(:,:,iS);
    else
        simX=W*SimFun+diagD.*simX+C+rtQ*simRand{iS};
    end
    sMeanTmp=sMeanTmp+simX;
    sCovTmp=sCovTmp+simX*simX';
    end
    simXC{SimLength+1}=simX;

%% REMARK: This version less prone to cancellation but slower
%simMean2=mean([simXC{2:end}],2);
%simCov2=([simXC{2:end}]-simMean)*([simXC{2:end}]-simMean)'/(SimLength*nSim-1);

%% REMARK: This version prone to cancellation when E[X] not approx. 0
simMean=sum(sMeanTmp,2)/(SimLength*nSim);
%% n-1 df for Cov, n df for Mean
simCov=(sCovTmp-(SimLength*nSim)*(simMean*simMean'))/(SimLength*nSim-1);

%disp([CorrR2_00(simMean,simMean2) CorrR2_00(simCov,simCov2)])


simCovOld=(1-decP)*simCov+decP*simCovOld;
Pbf=(1-decFix)*simCovOld+decFix*Pfix;
    else
%% 3.2 Otherwise Define Cov(X) based upon Prt estimates
Pbf=Prt*Prt'+Pfix;
    end
    end

for iBatch=1:nStart


%%%%%%%%%%%%%%


%% NOTE: WILL NEED TO CHANGE SetRep and iRepB for multiple data sets!!!!!!!!!!!!!!

%% Mapping between data-sets and the number of starts 
%% [will probably over-ride this later to equally-weight data as with MINDy]
SetRep=1:nRep;

ShiftBatch=0;

nStackEKF=numel(ShiftBatch);
do1X=nStackEKF==1;

%HrepRule
%dropEnd=kalSpec.nKal+kalSpec.nRec+1;
%[KdropSet,HdropSet]=BPK_Make_Hdrop(Hset,....;%% use BPK_gramschmidt..

%% 4 Minibatch specific setup

    iB=randi(size(MeasSet,2)-dropEnd);
    R=Rset;
%    Kdrop=KdropSet{Krule(iRepB)};
%    Hdrop=HdropSet{Hrule(iRepB)};
    ZX0=MeasSet(:,(iB));
    
    %% Optionally subtract influence of input on measurements
    if doInput
        ZX0=ZX0-BY*Uinput(:,(iB+ShiftBatch));
    end
    Hfull=Hset;
    Rtrue=Rset;
    Emat0=ErrMat;
    
    Imat=eye(size(MeasSet,1));

%Pbf=Prt*Prt'+Pfix;

%% 4.2 Initialize X0, P0 (note that this is ind. of measurement)
if ~strcmpi(KalSpec.BFcov,'n')

%    PP0=Pbf*Hfull';
%    tt0=PP0/(Hfull*PP0+Rtrue);

Sp=linsolve(Hfull*Pbf*Hfull'+Rtrue,Imat,opts);

BF=Pbf*Hfull'*Sp;
YSbf=ZX0'*Sp;
Xb0=BF*ZX0;

Gbf=(eye(nX)-BF*Hfull);
Pb0=Gbf*Pbf;
X=Xb0;
P=Pb0;
else
P=Pbase;
X=KalSpec.ScaleFun(iBatch)*Xest{iSet}(:,iB+ShiftBatch)+...
    (1-KalSpec.ScaleFun(iBatch))*Xguess{iSet}(:,iB+ShiftBatch);
end


%% 4.3 Select Measurement and Input data
Zkal=MeasSet{iSet}(:,iB+(1:(max(ShiftBatch)+(nStep))));
Zrec=MeasSet{iSet}(:,iB+(1:(max(ShiftBatch)+(nRec)))+nStep);

if doInput
    Ukal=Uinput{iSet}(:,iB+(0:(max(ShiftBatch)+(nStep))));
    Urec=Uinput{iSet}(:,iB+(0:(max(ShiftBatch)+(nRec)))+nStep);
    Zkal=Zkal-BY*Ukal(:,2:end);
    Zrec=Zrec-BY*Urec(:,2:end);
UK{nStep+1}=Ukal(:,nStep+ShiftBatch+1);
UR=cell(1,nRec+1);
for iRec=1:(nRec+1)
     UR{iRec}=Ukal(:,iRec+ShiftBatch);
end
end
H=Hset{iSet};

%% 4.4 Reset data cells
XC=cell(1,nStep);PC=cell(1,nStep);FunKal=cell(1,nStep);
JacC=cell(1,nStep);LC=cell(1,nStep);XrandPre=cell(1,nStep);Fchol=cell(1,nStep);
shiftXchol=cell(1,nStep);dEdFstep=cell(1,nStep);GC=cell(1,nStep);hYC=cell(1,nStep);
LpreGC=cell(1,nStep);Xrec=cell(1,nRec);FunRec=cell(1,nRec);dEdFrec=cell(1,nRec);
hYC2=cell(1,nStep);
d12x=1/(nX);
KT=cell(1,nStep);
Estore=zeros(1,nRec);
Ekal=zeros(1,nStep);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Forward Pass
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5 Kalman Filtering
for i1=1:nStep
    XC{i1}=X;
    PC{i1}=P;

    %% 5.1 Kalman Prediction
if isEKF
    Fun=tanh(SS.*X+V);
%    FunKal{i1}=Fun;
    X1t=W*Fun+diagD.*X+C;
    if do1X
        Fprime=SS.*(1-Fun.^2);
    else
    Fprime=SS.*sum(1-Fun.^2,2)/nStackEKF;
    end
    Jac=(W.*Fprime')+Dmat;
    P1t=Jac*P*Jac'+Q;
 %   JacC{i1}=Jac;
elseif isCKF
    LC=chol(P,'lower');
    RandPart=kron(LC,kSgn);
    Xp=(X+RandPart);
 %   XrandPre{i1}=Xp;
    Fun=-1+(2./(1+exp((-2*SS).*Xp-(2*V))));
%    Fchol{i1}=Fun;
    Xtmp=W*Fun+diagD.*Xp+C;
    X1t=Xtmp*BaseWeight';
    shiftXchol=Xtmp-X1t;
    P1t= shiftXchol*(shiftXchol)'/(2*nX)+Q;
elseif isUKF
    LC=chol(P,'lower');
    RandPart=LC*BaseMat;
    Xp=(X+RandPart);
  %  XrandPre{i1}=Xp;
    Fun=-1+(2./(1+exp((-2*SS).*Xp-(2*V))));
 %   Fchol{i1}=Fun;
    Xtmp=W*Fun+diagD.*(X+RandPart)+C;
    X1t=Xtmp*BaseWeight';
    shiftXchol=Xtmp-X1t;
    P1t=shiftXchol*(shiftXchol.*BaseWeight)'+Q;
end
    if doInput 
        X1t=X1t+Ukal(:,ShiftBatch+i1);
    end
%% 5.2 Measurement Prediction
    HP1C=H*P1t;
    yFull=Zkal(:,ShiftBatch+i1)-Hfull*X1t;
    tmp0=Emat0*yFull;
    Ekal(i1)=yFull(:)'*tmp0(:);
    y=Hproj*yFull;
%% 5.4 Kalman Gain
    tmp=(HP1C*H'+R)/2;
    S=linsolve(tmp+tmp',Imat,opts);
    K=HP1C'*S;
%% 5.5 Correction
    X=X1t+K*y;
   % if isNotEnd(i1)
   % GC{i1}=(I-K*H);    
    P=P1t-K*HP1C;
    P=(P+P')/2;
end

%% 6 Kalman-Free phase
%% 6.1 Predictions and Error
   for iRec=1:nRec
 %  Xrec{iRec}=X;
   Fun=-1+(2./(1+exp(-2*(SS.*X+V))));
   X=W*Fun+diagD.*X+C;
   if doInput
%   if doInputGradX
%       UR{iRec}=Ukal(:,iRec+ShiftBatch);
       X=X+BX*UR{iRec};
%   end
%   if doInputNOgradX
%       X=X+BX*Ukal(:,iRec+ShiftBatch);
   end
   tmpE=Zrec(:,iRec)-H*X;
   tmp0=Emat0*tmpE;
   Estore(iRec)=tmpE(:)'*tmp0(:);
   end
   Etotal=Etotal+Estore/NBatch;
end
end
