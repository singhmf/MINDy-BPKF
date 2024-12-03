function[Out]=BPKF_Full_tMask(Xguess,MeasSet,tMask,ParStr,KalSpec,ModelSpec,GradSpec,Uinput)
 startTIC=tic;
%% Remark: gradP isn't fully supported yet in terms of ParamStr and dependencies, Prt isn't updated?
%% tMask indicates which time-points to exclude from error function
VersStamp='3.26.24';
VersNotes='';

recSpace=0;

if isfield(ParStr,'PBT') && ~isempty(ParStr.PBT)
    ParStr.PBT=BPKF_Initialize_PBT(ParStr.PBT);
    PBTstr=ParStr.PBT;
    doPBT=true;
else
    doPBT=false;
end

[ParStr,BatchSchedule]=BPKF_Initialize_CV(ParStr);

Out.Settings={ParStr,KalSpec,ModelSpec,GradSpec};

%% 1: SETUP

%% 1.1 General Setup

%% 1.1A Ensure Dimensionality
if isfield(ModelSpec,'Wmask')
nX=size(ModelSpec.Wmask,1);
else
    if iscell(Xguess)
        nX=size(Xguess{1},1);
    else
        nX=size(Xguess,1);
    end
end
if ~iscell(MeasSet)
    MeasSet={MeasSet};
end
if ~iscell(Xguess)
    Xguess={Xguess};
end

if strcmpi(KalSpec.BFcov,'n')
Xest=Xguess;
end
%% 1.1.B Setup Identity matrix
I=eye(nX);
I2set=cell(size(MeasSet));
for ii=1:numel(MeasSet)
I2set{ii}=eye(size(MeasSet{ii},1));
end


%% 1.1C Determine if open-loop
if nargin<8 || isempty(Uinput)
    doInput=false;
else
    disp('Doing open-loop model')
    disp('REMEMBER TO UPDATE TO INCLUDE INFLUENCE OF U ON COV')
    doInput=true;
    if ~iscell(Uinput)
        Uinput={Uinput};
    end
    ModelSpec.nU=size(Uinput{1},1);
end
ModelSpec.doInput=doInput;

%% 1.2 Get Data Covariance (irrespective of input)
MeasCov=cell(size(MeasSet));
for ii=1:numel(MeasSet)
    MeasCov{ii}=cov(MeasSet{ii}');
end


%% 1.3 Model Initialization
[ParamValDep,ModelSpec,Hset,Rset,nadamStr]=BPKF_Initialize_Model(ModelSpec,GradSpec,numel(MeasSet));
ParamValDep.ParamStr.doInput=doInput;



probSet=@(yy)(cellfun(@(xx)(size(xx,2)-(KalSpec.nStep+KalSpec.nRec+1)),yy));

if isfield(ParStr,'CV')
markTrain=setdiff(1:numel(MeasSet),ParStr.CV.mark);
markCV=ParStr.CV.mark;
probSetCV=probSet(MeasSet(markCV));%/sum(probSet(markCV));
probSetCV=probSetCV/sum(probSetCV);
else
    markTrain=1:numel(MeasSet);
end
probSetTrain=probSet(MeasSet(markTrain));
probSetTrain=probSetTrain/sum(probSetTrain);

if ~isfield(ParStr,'CostFun')||isempty(ParStr.CostFun)
    ParStr.CostFun='L2';
    disp('Assuming L2 cost')
    costL2=true;
    doHuber=false;
else
    costL2=strcmpi(ParStr.CostFun,'L2');
    doHuber=strcmpi(ParStr.CostFun,'Huber');
end



if isfield(ParStr,'ParamValDep')
    disp('Reusing pre-specified start')
    ParamValDep=ParStr.ParamValDep;
end


if doInput
    doInputGradX=~(isfield(ParamValDep.ParamStr,'fixBX') && strcmpi(ParamValDep.ParamStr.fixBX(1),'y'));
    doInputNOgradX=~doInputGradX;
    doInputGradY=~(isfield(ParamValDep.ParamStr,'fixBY') && strcmpi(ParamValDep.ParamStr.fixBY(1),'y'));
else
    doInputGradX=false;
    doInputNOgradX=false;
    doInputGradY=false;
end

%% 1.4 More variable setup
opts.POSDEF = true;
opts.SYM = true;

BatchSz=ParStr.BatchSz;NBatch=ParStr.NBatch;

Pbase=KalSpec.Pbase;
%if isfield(ModelSpec,'nRep')
%    nRep=ModelSpec.nRep;
%else
%    nRep=1;
%end

minStep=KalSpec.minStep;
nStep=KalSpec.nStep;
nRec=KalSpec.nRec;
%% Scalefun is how quickly to interpolate between initial guess and past Kal estimates
%% for setting initial conditions on X
if ~isfield(KalSpec,'ScaleFun')
    disp('Assuming no ScaleFun')
    KalSpec.ScaleFun=@(xx)(0);
end
%% do Kal denotes when to skip doing the Kalman Filter
if ~isfield(KalSpec,'doKal')||mean(KalSpec.doKal==1,'all')==1
    doKal=true(1,nStep);
    isFixed_doKal=true;
elseif numel(KalSpec.doKal)==1
    isFixed_doKal=false;
    probKal=repmat(KalSpec.doKal,1,nStep);
elseif mean(ismember(KalSpec.doKal,[1 0]))==1&& numel(KalSpec.doKal)==nStep
    doKal=KalSpec.doKal;
    isFixed_doKal=true;
elseif numel(KalSpec.doKal)==nStep
    isFixed_doKal=false;
    probKal=KalSpec.doKal;
else
    error('Cannot interpret doKal--should be a scalar or vector of length nStep (or empty)')
end


%% 1.5 Setup error/cost matrix
ErrMat=cell(size(Hset));
if ~ischar(ModelSpec.ErrMat)
    ErrMat=BPKF_Convert_Cell_Size(ModelSpec.ErrMat,numel(Hset),'ErrMat');
elseif strcmpi(ModelSpec.ErrMat(1:2),'L2')
    for ii=1:numel(Hset)
        ErrMat{ii}=eye(size(Hset{ii},1));
    end
elseif strcmpi(ModelSpec.ErrMat(1:4),'Mahal')
    for ii=1:numel(Hset)
        ErrMat{ii}=inv(Hset{ii}*ParamStr.Qorig*Hset{ii}'+Rset{ii});
        ErrMat{ii}=(ErrMat{ii}+ErrMat{ii}')/2;
    end
end

%rtQ=ParamValDep.ParamStr.rtQ;
Prt=chol(Pbase,'lower');
%rtQ0=rtQ;

CovErr=ParStr.CovErr;

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
in1=1:nX;
in2=(nX+1):(2*nX);
TrilMat=tril(ones(nX))-eye(nX)/2;
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

%% 1.8 Setup recording parameters (for troubleshooting)
ParNames={'W','D','C','S','V','Q','BX','BY'};
recNames={};
for iP=1:numel(ParNames)
    tmpName=['rec',ParNames{iP}];
    if isfield(ParStr,tmpName)&&strcmpi(ParStr.(tmpName)(1),'y')
        recNames=[recNames,ParNames(iP)]; %#ok<AGROW> 
    end
end
nRec0=numel(recNames);
doSave=isfield(ParStr,'saveName');
if doSave
    doTmpSave=isfield(ParStr,'saveRate');
    if doTmpSave
        saveRate=ParStr.saveRate;
    end
else
    doTmpSave=false;
end

if doTmpSave
    save(ParStr.saveName,'Out')
%    matSAVE_Ahalf=matfile([ParStr.saveName,'_Ahalf'],'Writable',true);
%    matSAVE_Bhalf=matfile([ParStr.saveName,'_Bhalf'],'Writable',true);
end

if isfield(ParStr,'recPred')
    doRecPred=strcmpi(ParStr.recPred,'y');
else
    doRecPred=false;
end



if ~doInput
    recNames(strcmpi(recNames,'BX'))=[];
    recNames(strcmpi(recNames,'BY'))=[];
    if numel(recNames)~=nRec0
        disp('No input so not recording BX, BY')
    end
end

if isfield(ParStr,'recParam') && ~isempty(ParStr.recParam)
    recParam=ParStr.recParam;
else
    recParam=0;
end

if (~isfield(ParStr,'recSpace'))
    
    if numel(recNames)~=0
        disp('Assuming default spacing of 250 for ParStr.recSpace')
        recSpace=250;
    end
else
    recSpace=ParStr.recSpace;
end
if numel(recNames)==0
    recSpace=0;
end
if ~isempty(recNames)
    isQname=strcmpi(recNames,'Q');
end

if ~isfield(ParStr,'nStack')
    ShiftBatch=0;
else

%% Default spacing is 5 time-steps beween members of a stack minibatch
if isEKF && ~isempty(ParStr.nStack)
    if isfield(ParStr,'ShiftSpace')
        if isempty(ParStr.ShiftSpace)
            disp('ParStr.ShiftSpace is empty, using default of 5')
        else
        ShiftSpace=ParStr.ShiftSpace;
        end
    else
        ShiftSpace=5;
        disp('Using default spacing of 5 within EKF stacks')
    end
    ShiftBatch=ShiftSpace*((1:ParStr.nStack)-1);
else
    if ~isempty(ParStr.nStack)
    disp('Not doing EKF so ignoring nStack')
    end
    ShiftBatch=0;
end
end

Out.ShiftBatch=ShiftBatch;

if isfield(ModelSpec,'dropEnd')
    dropEnd=ModelSpec.dropEnd;
    if dropEnd<(nStep+1+nRec+max(ShiftBatch))
        disp('Specified dropEnd is too small for batch length--re-adjusting')
        dropEnd=nStep+1+nRec+max(ShiftBatch);
    end
else
    dropEnd=nStep+1+nRec+max(ShiftBatch);
end
nStackEKF=numel(ShiftBatch);
do1X=nStackEKF==1;

doFullDistr=isfield(ParStr,'saveDistr')&&~isempty(ParStr.saveDistr);
if doFullDistr
    saveDistrRate=ParStr.saveDistr;
if do1X
    if isfield(ParStr,'CV')&&~isempty(ParStr.CV)
    distrK={zeros(nStep,nStackEKF),zeros(nStep,ParStr.CV.Size)};
    distrR={zeros(nRec,ParStr.BatchSz),zeros(nRec,ParStr.CV.Size)};
    else
        distrK={zeros(nStep,nStackEKF)};
        distrR={zeros(nRec,nStackEKF)};
    end
else
       disp('Using batch KF so saving single-batch distributions')
       distrK=repmat({zeros(nStep,nStackEKF)},1,2);
       distrR=repmat({zeros(nRec,nStackEKF)},1,2);
end 

     Out.DistrKal=cell(1,ceil(NBatch/saveDistrRate));
     Out.DistrRec=cell(1,ceil(NBatch/saveDistrRate));
     if isfield(ParStr,'CV')
         Out.cvDistrKal=cell(1,abs(min(BatchSchedule)));
         Out.cvDistrRec=cell(1,abs(min(BatchSchedule)));
     end
end

%% Total evaluations (denominator in mean-error)
SzScale=BatchSz*(nRec+nStep)*numel(ShiftBatch);


doAdaR=isfield(KalSpec,'alphaR');
doAdaQ=isfield(KalSpec,'alphaQ');
%% Adaptive KF Setup
if doAdaR
    posR=repmat(reshape(Rset,numel(Rset),1),1,nStep);
    Rbase=Rset;
    alphaR=KalSpec.alphaR;
    betaR=KalSpec.betaR;
    mixR=KalSpec.mixR;
    mixR(3)=1-(mixR(2)+mixR(1));
end
if doAdaQ
    posQ=repmat({ModelSpec.Q},1,nStep);
    alphaQ=KalSpec.alphaQ;
    betaQ=KalSpec.betaQ;
    mixQ=KalSpec.mixQ;
    mixQ(3)=1-(mixQ(1)+mixQ(2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Out.StartupTime=toc(startTIC);

Out.BatchSchedule=BatchSchedule;

if doHuber
    decHuber=ParStr.decHuber;
    HuberScale=ParStr.HuberScale;
    ParamValDep.decHuber=decHuber;
    ParamValDep.HuberScale=HuberScale;
  %  Out.kalE_Huber=Out.kalE;
  %  Out.recE_Huber=Out.recE;
  %  if isfield(Out,'recCV')
  %      Out.kalCV_Huber=Out.kalCV;
  %      Out.recCV_Huber=Out.recCV;
  %  end

    med_recHuber=zeros(nRec,ParStr.BatchSz);
    med_kalHuber=zeros(nStep,ParStr.BatchSz);
    E_kalHuber0=zeros(nStep,1);
    E_recHuber0=zeros(nRec,1);

    med_Huber_recSet=cell(size(MeasSet));
    med_Huber_kalSet=cell(size(MeasSet));
    Etot=zeros(size(MeasSet));
    for iH=1:numel(Etot)
        Etot(iH)=single(HuberScale*median(sum(MeasSet{iH}.*(ErrMat{iH}*MeasSet{iH}),1)));
        med_Huber_recSet{iH}=repmat(Etot(iH),nRec,1);
        med_Huber_kalSet{iH}=repmat(Etot(iH),nStep,1);
    end

end
simXC=cell(1,SimLength+1);


if strcmpi(KalSpec.BFcov(1),'s')&& CovErr~=0
    Out.CovErr=zeros(1,max(BatchSchedule));
end

Out.recE=nan(nRec,max(BatchSchedule));
Out.kalE=nan(nStep,max(BatchSchedule));

Out.recE_Huber=nan(nRec,max(BatchSchedule));
Out.kalE_Huber=nan(nStep,max(BatchSchedule));


Out.recCV=nan(nRec,abs(min(BatchSchedule)));
Out.kalCV=nan(nStep,abs(min(BatchSchedule)));

Out.recCV_Huber=nan(nRec,abs(min(BatchSchedule)));
Out.kalCV_Huber=nan(nStep,abs(min(BatchSchedule)));

isLIN=isfield(ModelSpec,'SigTy')&& strcmpi(ModelSpec.SigTy(1:3),'Lin');
if isLIN
    disp('Doing LINEAR: Ignore S, V')
    isEKF=false;
    isCKF=false;
    isUKF=false;
end

batchTIC=tic;
for iBatch=BatchSchedule
%disp(iBatch)
    %% 2 Initialization
if doHuber
    E_kalHuber0=zeros(nStep,1);
    E_recHuber0=zeros(nRec,1);
end
    %% 2.1 Clear gradients
    gradW=0;gradD=0;gradSS=0;gradC=0;gradV=0;gradPrt=0;gradQrt=0;
    gradBY=0;gradBX=0;
    dEdPbf=0;
    UK=cell(1,nStep+1);
    Estore0=zeros(1,nRec);
    Ekal0=zeros(1,nStep);

    Hproj=1;
    if iBatch>0
        isCVrun=false;
    if mod(iBatch,250)==1
        disp([iBatch NBatch]);
        Out.FinBatch=iBatch;
        Out.runTime=toc(batchTIC);
    end
    else
        isCVrun=true;
    end
        
    %% 2.2 Apply Dependencies
    if doInput
    [W,diagD,C,SS,V,rtQ,BX,BY]=BPKF_ExtractDependency(ParamValDep);
    else
    [W,diagD,C,SS,V,rtQ]=BPKF_ExtractDependency(ParamValDep);
    end
    Dmat=diag(diagD);

    ATS=SS.*(W');

    nX2=2*nX;
    n2V=-2*V;    
    n2S=-2*SS;
%% BF Initialization
%    if ~strcmpi(ModelSpec.freeQ(1),'n')
%    rtQ=rtQ0+ModelSpec.baseQchol;
isDiagQ=min(size(rtQ))==1;
if isDiagQ
    dgRtQ=rtQ;
    rtQ=diag(rtQ);
    Q=rtQ.^2;
else
    Q=rtQ*rtQ';
end
%    end

%% 3 Initialize distributions
    if ~strcmpi(KalSpec.BFcov(1),'n')
        %% 3.1 Use simulations (if specified)
    if strcmpi(KalSpec.BFcov(1),'s')
       
    simX=simStart(:,randperm(nSaveStart,nSim));
%    simRand=randn(nX,nSim,SimLength);
    simFunC=cell(1,SimLength);
    sMeanTmp=0*sMeanTmp;
    sCovTmp=0*sCovTmp;
    for iS=1:SimLength
    simRand{iS}=randn(nX,nSim);
    simXC{iS}=simX;
%    SimFun=-1+2./(1+exp(-2*(SS.*simX+V)));
if isLIN
    SimFun=simX;
else
    SimFun=-1+2./(1+exp(n2S.*simX+n2V));
end
    simFunC{iS}=SimFun;
    if isDiagQ
        simX=W*SimFun+diagD.*simX+C+dgRtQ.*simRand{iS};%(:,:,iS);
    else
        simX=W*SimFun+diagD.*simX+C+rtQ*simRand{iS};
    end
 %   sMeanTmp=sMeanTmp+simX;
 %   sCovTmp=sCovTmp+simX*simX';
    end
    simXC{SimLength+1}=simX;

%% REMARK: This version less prone to cancellation but slower
simMean=mean([simXC{2:end}],2);
simCov=([simXC{2:end}]-simMean)*([simXC{2:end}]-simMean)'/(SimLength*nSim-1);

%% REMARK: This version prone to cancellation when E[X] not approx. 0
%simMean=sum(sMeanTmp,2)/(SimLength*nSim);
%% n-1 df for Cov, n df for Mean
%simCov=(sCovTmp-(SimLength*nSim)*(simMean*simMean'))/(SimLength*nSim-1);

%disp([CorrR2_00(simMean,simMean2) CorrR2_00(simCov,simCov2)])


simCovOld=(1-decP)*simCov+decP*simCovOld;
Pbf=(1-decFix)*simCovOld+decFix*Pfix;
    else
%% 3.2 Otherwise Define Cov(X) based upon Prt estimates
Pbf=Prt*Prt'+Pfix;
    end
    end
%%%%%%%%%%%%%%


%% Mapping between data-sets and the number of starts 
if isCVrun
SetRep=BPKF_Discrete_Sample(ParStr.CV.Size,probSetCV);
SetRep=markCV(SetRep);
else
SetRep=BPKF_Discrete_Sample(ParStr.BatchSz,probSetTrain);
SetRep=markTrain(SetRep);
end


%HrepRule
%dropEnd=kalSpec.nKal+kalSpec.nRec+1;
%[KdropSet,HdropSet]=BPK_Make_Hdrop(Hset,....;%% use BPK_gramschmidt..


isDistr0=doFullDistr&&(or(isCVrun,mod(iBatch,saveDistrRate)==1));

if isDistr0
    kalDistr=distrK{1+isCVrun};
    recDistr=distrR{1+isCVrun};
end



%% 4 sample specific setup
for iMiniBatch=1:numel(SetRep)

    if isDistr0
        if do1X
            isDistrBatch=true;
        else
            isDistrBatch=iMiniBatch==1;
        end
    else
        isDistrBatch=false;
    end



    %% 4.1 Measurement Specification
    iSet=SetRep(iMiniBatch);
    iB=randi(size(MeasSet{iSet},2)-dropEnd);
    if doAdaR
    R0=Rset{iSet};
    tmpR=Rbase{iSet};
    else
    R=Rset{iSet};
    end
    if doAdaQ
    Q0=Q;
    tmpQ=posQ{1};
    end


    %% randomly assign which steps use the Kalman Filter (optional)
    if ~isFixed_doKal
        doKal=rand(1,nStep)<=probKal;
    end
    
%    Kdrop=KdropSet{Krule(iRepB)};
%    Hdrop=HdropSet{Hrule(iRepB)};
    ZX0=MeasSet{iSet}(:,(iB+ShiftBatch));
    
    %% Optionally subtract influence of input on measurements
    if doInput
        ZX0=(ZX0-BY*Uinput{iSet}(:,(iB+ShiftBatch))).*tMask{iSet}(:,(iB+ShiftBatch));
    end
    Hfull=Hset{iSet};
    Rtrue=Rset{iSet};
    Emat0=ErrMat{iSet};
 
    HC20=-2*Hfull'*ErrMat{iSet};
    
    if doInput
        HCinput=-2*ErrMat{iSet};
    end
    Imat=I2set{iSet};

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

tMaskKal=tMask{iSet}(:,iB+(1:(max(ShiftBatch)+(nStep))));
tMaskRec=tMask{iSet}(:,iB+(1:(max(ShiftBatch)+(nRec)))+nStep);

if doInput
    Ukal=Uinput{iSet}(:,iB+(0:(max(ShiftBatch)+(nStep))));
    Urec=Uinput{iSet}(:,iB+(0:(max(ShiftBatch)+(nRec)))+nStep);
    Zkal=Zkal-BY*Ukal(:,2:end);
    Zrec=Zrec-BY*Urec(:,2:end);
UK{nStep+1}=Ukal(:,nStep+ShiftBatch+1);
UR=cell(1,nRec+1);
for iRec=1:(nRec+1)
     UR{iRec}=Urec(:,iRec+ShiftBatch);
end
end
H=Hset{iSet};
if doHuber
kalHuberC=med_Huber_kalSet{iSet};
recHuberC=med_Huber_recSet{iSet};
end
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

E_recHuber=zeros(nRec,1);
E_kalHuber=zeros(nStep,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Forward Pass
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5 Kalman Filtering
for iKal=1:nStep
    XC{iKal}=X;
    PC{iKal}=P;
if doAdaR
    R=mixR(1)*tmpR+mixR(2)*posR{iSet,iKal}+mixR(3)*R0;
    R=(R+R')/2;
end
if doAdaQ
    Q=mixQ(1)*tmpQ+mixQ(2)*posQ{iKal}+mixQ(3)*Q0;
    Q=(Q+Q')/2;
end
    %% 5.1 Kalman Prediction
if isEKF
    Fun=-1+(2./(1+exp((-2*SS).*X-(2*V))));
    
    FunKal{iKal}=Fun;
    X1t=W*Fun+diagD.*X+C;
    if do1X
        Fprime=SS.*(1-Fun.^2);
    else
    Fprime=SS.*sum(1-Fun.^2,2)/nStackEKF;
    end
    Jac=(W.*Fprime')+Dmat;
    P1t=Jac*P*Jac'+Q;
    JacC{iKal}=Jac;
elseif isCKF
    LC{iKal}=chol(P,'lower');
    RandPart=kron(LC{iKal},kSgn);
    Xp=(X+RandPart);
    XrandPre{iKal}=Xp;
    Fun=-1+(2./(1+exp((-2*SS).*Xp-(2*V))));
    Fchol{iKal}=Fun;
    Xtmp=W*Fun+diagD.*Xp+C;
    X1t=Xtmp*BaseWeight';
    shiftXchol{iKal}=Xtmp-X1t;
    P1t= shiftXchol{iKal}*(shiftXchol{iKal})'/(2*nX)+Q;
elseif isUKF
    LC{iKal}=chol(P,'lower');
    RandPart=LC{iKal}*BaseMat;
    Xp=(X+RandPart);
    XrandPre{iKal}=Xp;
    Fun=-1+(2./(1+exp((-2*SS).*Xp-(2*V))));
    Fchol{iKal}=Fun;
    Xtmp=W*Fun+diagD.*(X+RandPart)+C;
    X1t=Xtmp*BaseWeight';
    shiftXchol{iKal}=Xtmp-X1t;
    P1t=shiftXchol{iKal}*(shiftXchol{iKal}.*BaseWeight)'+Q;
elseif isLIN
    Fun=X;
    FunKal{iKal}=Fun;
    X1t=W*Fun+diagD.*X+C;
    Jac=W+Dmat;
    P1t=Jac*P*Jac'+Q;
    JacC{iKal}=Jac;

end
    if doInputGradX
        UK{iKal}=Ukal(:,ShiftBatch+iKal);
        X1t=X1t+BX*UK{iKal};
    end
    if doInputNOgradX %% Storing Inputs in UC is only needed for gradients
        X1t=X1t+BX*Ukal(:,ShiftBatch+iKal);
    end
%% 5.2 Measurement Prediction
    HP1C=H*P1t;

    tMask0=tMaskKal(iKal+ShiftBatch);
    tMean0=mean(tMask0);
    yFull=(Zkal(:,ShiftBatch+iKal)-Hfull*X1t).*tMask0;
    tmp0=Emat0*yFull/((tMean0==0)+tMean0);
    Ekal(iKal)=yFull(:)'*tmp0(:);
    y=Hproj*yFull;
%% 5.3 Error calculation
 %   if i1
    dEdFstep{iKal}=(HC20*yFull)/((tMean0==0)+tMean0);
   if doHuber
       HuberTmp=sqrt(1+(sum(yFull.*tmp0,1)/kalHuberC(iKal)));%/cHuber);
       E_kalHuber(iKal)=2*kalHuberC(iKal)*(-1+mean(HuberTmp));
       dEdFstep{iKal}=dEdFstep{iKal}./HuberTmp;
       if ~isCVrun
       med_kalHuber(iKal,iMiniBatch)=median(HuberTmp);
       end
   end

    if isDistrBatch

        if do1X
            kalDistr(iKal,iMiniBatch)=(yFull)'*(Emat0*yFull);
        else
            kalDistr(iKal,:)=sum(yFull.*(Emat0*yFull),1);
        end
    end
 %   end
    if doInputGradY
        if iKal>minStep
            gradBY=gradBY+HCinput*yFull*(Ukal(:,ShiftBatch+iKal+1))';
        end
    end
    if doKal(iKal)
%% 5.4 Kalman Gain
%    disp([iBatch iKal])
    halfS=(HP1C*H'+R)/2;
    %try
    invS=linsolve(halfS+halfS',Imat,opts);
    %catch
   %     disp([iKal iKal])
   %     disp(eig(halfS+halfS'))
    %end
    hYC{iKal}=H'*(invS*y);  
    hYC2{iKal}=hYC{iKal}/2;
    K=HP1C'*invS;
if doInputGradY
        KT{iKal}=K';
end
%% 5.5 Correction
    X=X1t+K*y;
   % if isNotEnd(i1)
    GC{iKal}=(I-K*H);    
    P=P1t-K*HP1C;
    if doAdaR
       deltaR=(yFull-H*(X-X1t))*(yFull-H*(X-X1t))'/nStackEKF+HP1C*H';
       posR{iSet,iKal}=(1-alphaR)*posR{iSet,iKal}+(alphaR)*deltaR;
       tmpR=(1-betaR)*tmpR+betaR*deltaR;
    end
    if doAdaQ
       deltaQ=(X-X1t)*(X-X1t)'/nStackEKF;
       posQ{iKal}=(1-alphaQ)*posQ{iKal}+(alphaQ)*deltaQ;
       tmpQ=(1-betaQ)*tmpQ+betaQ*deltaQ;
    end



    else
        P=P1t;
        X=X1t;
    end
    P=(P+P')/2;
    if isCKF && (iKal~=1)&&(doKal(iKal-1))
    LpreGC{iKal-1}=LC{iKal}\GC{iKal-1};
    end
end
if strcmpi(KalSpec.BFcov(1),'s')
simStart(:,randi(nSaveStart))=X(:,1);
end
if doAdaR
    Rbase{iSet}=.2*tmpR+.8*Rbase{iSet};
end
%% 6 Kalman-Free phase
%% 6.1 Predictions and Error
   for iRec=1:nRec
   Xrec{iRec}=X;
   if isLIN
       Fun=X;
   else
   Fun=-1+(2./(1+exp(-2*(SS.*X+V))));
   end
   X=W*Fun+diagD.*X+C;
   FunRec{iRec}=Fun;
   if doInput
%   if doInputGradX
%       UR{iRec}=Ukal(:,iRec+ShiftBatch);
       X=X+BX*UR{iRec};
%   end
%   if doInputNOgradX
%       X=X+BX*Ukal(:,iRec+ShiftBatch);
   end
   
   tMask0=(tMaskRec(iRec+ShiftBatch));
   tMean0=mean(tMask0);
   tmpE=(Zrec(:,iRec+ShiftBatch)-H*X).*tMask0;
   dEdFrec{iRec}=HC20*tmpE/((tMean0==0)+tMean0);
   
   tmp0=Emat0*tmpE;
   Estore(iRec)=tmpE(:)'*tmp0(:)/((tMean0==0)+tMean0);
   if doHuber
       HuberTmp=sqrt(1+sum(tmpE.*tmp0,1)/recHuberC(iRec));%/cHuber);
       E_recHuber(iRec)=2*recHuberC(iRec)*(-1+mean(HuberTmp));
       dEdFrec{iRec}=dEdFrec{iRec}./HuberTmp;
       if ~isCVrun
       med_recHuber(iRec,iMiniBatch)=median(HuberTmp);
       end
   end
   if isDistrBatch
       if do1X
           recDistr(iRec,iMiniBatch)=tmpE'*tmp0;
       else
           recDistr(iRec,:)=sum(tmpE.*tmp0,1);
       end
   end
   if doInputGradY
       gradBY=gradBY+HCinput*tmpE*UR{iRec+1}';
   end
   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Back Prop Phase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isCVrun
%% 6.2 Free-sim Backprop
    dEdF=0;
    for iRec=nRec:-1:1
    dEdF=dEdF+dEdFrec{iRec};
    X=Xrec{iRec};
    Fun=FunRec{iRec};
gradW=gradW+dEdF*Fun';
gradD=gradD+sum(dEdF.*X,2);
if doInputGradX
gradBX=gradBX+dEdF*UR{iRec}';
end
gradC=gradC+sum(dEdF,2);
if isLIN
    dEdTan=(W'*dEdF);
dEdF=dEdTan+dEdF.*diagD;
else
    dEdTan=(W'*dEdF).*(1-Fun.^2);
gradSS=gradSS+sum(dEdTan.*X,2);
gradV=gradV+sum(dEdTan,2);
dEdF=SS.*dEdTan+dEdF.*diagD;
end
    end
dEdX=dEdF;
dEdP=zeros(nX);

%% 7. BP EKF

if isEKF
    %% 7.1: BP to predicted State
    for ii=nStep:-1:1
        if doKal(ii)
    GcEx=GC{ii}'*dEdX;  
        else
            GcEx=dEdX;
        end
    if(ii>minStep)
    GcEx2=GcEx+dEdFstep{ii};
    else
        GcEx2=GcEx;
    end

    %% 7.2 Parameter Gradients w.r.t. state

Fun=FunKal{ii};

dTan=(1-Fun.^2);
    %% REMARK: This assumes nStackEKF=1
dEdTan=(W'*GcEx2).*dTan;

gradW=gradW+GcEx2*Fun';
gradD=gradD+sum(XC{ii}.*GcEx2,2);
gradC=gradC+sum(GcEx2,2);
gradSS=gradSS+sum(XC{ii}.*dEdTan,2);
gradV=gradV+sum(dEdTan,2);


%% 7.3 Backprop to predicted covariance
%dEdPt1=hYC{ii}*GcEx'+(GC{ii}'*dEdP*GC{ii});
%gradQrt=gradQrt+(dEdPt1+dEdPt1')*rtQ;
%dEdPt1Jac2=(dEdPt1+dEdPt1')*JacC{ii};
if doKal(ii)
    if ii==nStep
       dEdPt1=hYC{ii}*GcEx'+GcEx*hYC{ii}';
    else
        dEdPt1=hYC{ii}*GcEx'+GcEx*hYC{ii}'+2*(GC{ii}'*dEdP*GC{ii});
    end
else
dEdPt1=dEdP;
end
gradQrt=gradQrt+(dEdPt1)*rtQ;
dEdPt1Jac2=(dEdPt1)*JacC{ii};


tmpPGF=dEdPt1Jac2*(PC{ii});
if doInputGradX
    gradBX=gradBX+GcEx2*UK{ii}';
end
if doInputGradY && doKal(ii)
    gradBY=gradBY-(KT{ii}*dEdX)*UK{ii+1}';
end
%% 7.4 Parameter gradients with respect to covariance
Fprime=SS.*dTan;
if do1X
gradW=gradW+tmpPGF.*(Fprime');
else
gradW=gradW+(tmpPGF).*mean(Fprime,2)';    
end
gradD=gradD+diag(tmpPGF);

dJacdTan=sum((W.*tmpPGF),1)';
dF2dXV=(-2*Fun.*Fprime)/nStackEKF;
dF2dXdS=(1-2*SS.*XC{ii}.*Fun).*dTan/nStackEKF;
gradSS=gradSS+sum(dJacdTan.*dF2dXdS,2);
gradV=gradV+sum(dJacdTan.*dF2dXV,2);
%if ii~=1
dEdP=JacC{ii}'*dEdPt1Jac2/2;
dEdX=SS.*dJacdTan.*dF2dXV+SS.*dEdTan+diagD.*GcEx2;
%end
    end
%% 8. BP-CKF (instead of 7)
elseif isCKF
    dEdPpreL=zeros(nX);
for ii=nStep:-1:1   
    
    %% 8.1 BP to prior distribution
    if doKal(ii)
    GcEx=GC{ii}'*dEdX;    
    else
        GcEx=dEdX;
    end
    if(ii>minStep)
    GcEx2=GcEx+dEdFstep{ii};
    else
        GcEx2=GcEx;
    end
  %  Fun=FunKal{ii};%tanh(SS.*XC{ii});
        dFun=1-Fchol{ii}.^2;
%        dEdPt1=hYC{ii}*GcEx'+(LpreGC{ii}'*dEdPpreL*LpreGC{ii});
%        dEPdX=(dEdPt1+dEdPt1')*(shiftXchol{ii}.*BaseWeight);
if doKal(ii)
if ii==nStep
        dEdPt1=GcEx*hYC2{ii}'+hYC2{ii}*GcEx';
else
        dEdPt1=GcEx*hYC2{ii}'+hYC2{ii}*GcEx'+LpreGC{ii}'*dEdPpreL*LpreGC{ii};
end
else
    if ii==nStep
        dEdPt1=dEdP;
    else
        invLC=linsolve(LC{ii+1},InX,opts2);
        dEdPt1=invLC'*dEdPpreL*invLC; %% =dEdP  
    end
end
        gradQrt=gradQrt+dEdPt1*rtQ;
        %% Remark: This line is only true for CKF (not UKF)
        dEPdX=dEdPt1*(shiftXchol{ii}.*d12x);
 %               dEPdX=dEdPt2*(shiftXchol{ii}*d12x);
%        dEPdX2=(dEdPt1+dEdPt1')*(shiftXchol{ii})/(2*nX);       
     %   dEdXchol=dEPdX+(GcEx2).*BaseWeight;
        dEdXchol=dEPdX+(GcEx2./nX2);
        %% 8.2 BP to Parameter
        if doInputGradX
            gradBX=gradBX+dEdX;
        end
       
        gradC=gradC+sum(dEdXchol,2);
        gradD=gradD+sum(dEdXchol.*XC{ii},2);
        gradW=gradW+dEdXchol*Fchol{ii}';
        dEdXpre=(dFun).*(ATS*dEdXchol)+diagD.*dEdXchol;
%                dEdXpre=(((SS.*dFun).))*dEdXchol;%+diagD.*dEdXchol;        
        dEdX=sum(dEdXpre,2);
     %   dB=dEdXpre*BaseMat';
%        dB=dEdXpre*BaseMatT;
        dB=dEdXpre(:,in1)-dEdXpre(:,in2);
    QZ=LC{ii}'*dB;
%    invLC=pinv(LC{ii});
    dEdPpreL=TrilMat.*QZ;%(tril(QZ)-diag(diag(QZ))/2);
    dEdPpreL=(dEdPpreL+dEdPpreL')/2;
%    tmp=invLC'*(t0)*invLC;
%    tmp=(LC{ii}'\t0)/LC{ii};
%    dEdP=(tmp+tmp')/2;
end
    %% 8.3 BP to P0
    invLC=linsolve(LC{ii},InX,opts2);
%    invLC=pinv(LC{ii});
   t0=TrilMat.*QZ;%(tril(QZ)-diag(diag(QZ))/2);
    tmp=invLC'*(t0)*invLC;
 %   tmp=(LC{ii}'\t0)/LC{ii};
    dEdP=(tmp+tmp')/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif isLIN
    %% 7.1: BP to predicted State
    for ii=nStep:-1:1
        if doKal(ii)
    GcEx=GC{ii}'*dEdX;  
        else
            GcEx=dEdX;
        end
    if(ii>minStep)
    GcEx2=GcEx+dEdFstep{ii};
    else
        GcEx2=GcEx;
    end

    %% 7.2 Parameter Gradients w.r.t. state

Fun=FunKal{ii};

    %% REMARK: This assumes nStackEKF=1
dEdTan=(W'*GcEx2);

gradW=gradW+GcEx2*Fun';
gradD=gradD+sum(XC{ii}.*GcEx2,2);
gradC=gradC+sum(GcEx2,2);


%% 7.3 Backprop to predicted covariance
%dEdPt1=hYC{ii}*GcEx'+(GC{ii}'*dEdP*GC{ii});
%gradQrt=gradQrt+(dEdPt1+dEdPt1')*rtQ;
%dEdPt1Jac2=(dEdPt1+dEdPt1')*JacC{ii};
if doKal(ii)
    if ii==nStep
       dEdPt1=hYC{ii}*GcEx'+GcEx*hYC{ii}';
    else
        dEdPt1=hYC{ii}*GcEx'+GcEx*hYC{ii}'+2*(GC{ii}'*dEdP*GC{ii});
    end
else
dEdPt1=dEdP;
end
gradQrt=gradQrt+(dEdPt1)*rtQ;
dEdPt1Jac2=(dEdPt1)*JacC{ii};


tmpPGF=dEdPt1Jac2*(PC{ii});
if doInputGradX
    gradBX=gradBX+GcEx2*UK{ii}';
end
if doInputGradY && doKal(ii)
    gradBY=gradBY-(KT{ii}*dEdX)*UK{ii+1}';
end
%% 7.4 Parameter gradients with respect to covariance

gradW=gradW+(tmpPGF);
gradD=gradD+diag(tmpPGF);

%if ii~=1
dEdP=JacC{ii}'*dEdPt1Jac2/2;
dEdX=dEdTan+diagD.*GcEx2;
%end
    end
end

%% 9. Backprop over global distributions
if ~strcmpi(KalSpec.BFcov(1),'n')

%% 9.1 BP over least-squares
if doInputGradY
    gradBY=gradBY-(BF'*dEdX)*UK{1}';
end

tmp=H'*(YSbf')*(Gbf'*dEdX)';
dEdPbf=dEdPbf+(tmp+tmp')/2+Gbf'*dEdP*Gbf;
gradPrt=gradPrt+2*dEdPbf*Prt; 
if strcmpi(KalSpec.BFcov(1),'s')&& CovErr~=0
    if iMiniBatch==1
Out.CovErr(iBatch)=sum((MeasCov{iSet}-H*simCovOld*H'-R).^2,[1 2]);
    end
dEdPbf=dEdPbf-(2*CovErr/(1-decFix))*H'*(MeasCov{iSet}-H*simCovOld*H'-R)*H;
end
end

end
%% Accumulate error between batches
Estore0=Estore0+Estore;
Ekal0=Ekal+Ekal0;

if doHuber
    E_kalHuber0=E_kalHuber0+E_kalHuber;
    E_recHuber0=E_recHuber0+E_recHuber;
end

%% end of minibatch sample
end



%% 9.2 BP over simulation
if strcmpi(KalSpec.BFcov(1),'s')
dEdF=zeros(nX,nSim);
zScale=(2/(SimLength*nSim-1))*(1-decP)*(1-decFix);
    for iS=SimLength:-1:1 
            dEdF=dEdF+dEdPbf*(zScale*(simXC{iS+1}-simMean)); 
          %  dEdF=dEdF+zScale*dEdPbf*(simXC{iS+1}-simMean);  
   X=simXC{iS};


gradQrt=gradQrt+dEdF*simRand{iS}';
%gradQrt=gradQrt+dEdF*simRand(:,:,iS)';
 %   Fun=tanh(SS.*X+V);%FunRec{iRec};
    Fun=simFunC{iS};
gradW=gradW+dEdF*Fun';
 gradD=gradD+sum(dEdF.*X,2);
 gradC=gradC+sum(dEdF,2);
 if isLIN
dEdTan=(W'*dEdF);
 dEdF=((dEdTan))+diagD.*dEdF;
 else
dEdTan=(W'*dEdF).*(1-Fun.^2);
 gradV=gradV+sum(dEdTan,2);
 gradSS=gradSS+sum(dEdTan.*X,2);
 dEdF=SS.*((dEdTan))+diagD.*dEdF;
 end
    end
end



%% 10. Record Diagnostics
%% 10.1 Store mean error
if isCVrun
Out.recCV(:,abs(iBatch))=Estore0/(ParStr.CV.Size*nStackEKF);
Out.kalCV(:,abs(iBatch))=Ekal0/(ParStr.CV.Size*nStackEKF);

if doHuber
    Out.recCV_Huber(:,abs(iBatch))=E_recHuber0/(ParStr.CV.Size*nStackEKF);
    Out.kalCV_Huber(:,abs(iBatch))=E_kalHuber0/(ParStr.CV.Size*nStackEKF);
end

Out.batchCV(:,abs(iBatch))=BatchSchedule(max(1,-1+find(BatchSchedule==iBatch)));

if doFullDistr
    Out.cvDistrKal{abs(iBatch)}=kalDistr;
    Out.cvDistrRec{abs(iBatch)}=recDistr;
    Out.cvInd{1,abs(iBatch)}=[SetRep, iB];
end


else

if doFullDistr&&mod(iBatch,saveDistrRate)==1
    Out.DistrKal{ceil(iBatch/saveDistrRate)}=kalDistr;
    Out.DistrRec{ceil(iBatch/saveDistrRate)}=recDistr;
    Out.DistrInd{1,ceil(iBatch/saveDistrRate)}=[SetRep, iB];
end


Out.recE(:,iBatch)=Estore0/(nStackEKF*BatchSz);
Out.kalE(:,iBatch)=Ekal0/(nStackEKF*BatchSz);


if doHuber
Out.recE_Huber(:,iBatch)=E_recHuber0/(nStackEKF*BatchSz);
Out.kalE_Huber(:,iBatch)=E_kalHuber0/(nStackEKF*BatchSz);
end

%% 10.3 Record parameter estimates (if specified)
%% REMARK: Recorded estimates are those used on this batch (not after update)

%%%%%%%%%%%%%%%%%%%%%
%% Record Parameters
%%%%%%%%%%%%%%%%%%%%%
if (recParam~=0) && mod(iBatch,recParam)==1
    recordInd=ceil(iBatch/recParam);
 %   ParTmp=BPKF_Apply_Dependency(ParamValDep);
        Out.recParam{recordInd}=ParamValDep;
        Out.recParamInd(1,recordInd)=iBatch;
end
if (recSpace~=0) && mod(iBatch,recSpace)==1
    recordInd=ceil(iBatch/recSpace);
    ParTmp=BPKF_Apply_Dependency(ParamValDep);
    if iBatch==1
    for iP=1:numel(recNames)
        tmpName=strcat('rec',recNames{iP});
        if isQname(iP)
            Out.(tmpName)=zeros(numel(ParTmp.Qchol),ceil(NBatch/recSpace));
        else
            Out.(tmpName)=zeros(numel(ParTmp.(recNames{iP})),ceil(NBatch/recSpace));
        end
    end
    end
    disp('recording')
    for iP=1:numel(recNames)
        tmpName=strcat('rec',recNames{iP});
        if isQname(iP)
            Out.(tmpName)(:,recordInd)=ParTmp.Qchol(:);
        else
            Out.(tmpName)(:,recordInd)=ParTmp.(recNames{iP})(:);
        end
    end
 %   if doRecSD
 %       Out.recSD(:,recordInd)=median(recSD,2);
 %       Out.kalSD(:,recordInd)=median(kalSD,2);
 %   end

end
if doTmpSave && mod(iBatch,saveRate)==1 && iBatch>1
    
[Out.Param,Out.rtQ]=BPKF_Condense_Param(ParamValDep,doInput);

if strcmpi(KalSpec.BFcov(1),'s')
    Out.Cov=simCovOld;
    Out.InstCov=simCov;
    Out.MeasCov=MeasCov;
end
disp('Saving')
ttSave=tic;
%if mod(iBatch,saveRate*2)==1
%matSAVE_Bhalf.Out=Out;
%else
%matSAVE_Ahalf.Out=Out;
%end
save(ParStr.saveName,'Out')
toc(ttSave)
disp(['Done saving: ',num2str(toc(ttSave))])
end


%% 11. Parameter updates

%% 11.1 Rescale gradients based upon batchsz and negate
GradVal.W=-gradW/SzScale;
GradVal.D=-gradD/SzScale;
GradVal.C=-gradC/SzScale;
GradVal.S=-gradSS/SzScale;
GradVal.V=-gradV/SzScale;
if isLIN
    GradVal.S=0*GradVal.S;
    GradVal.V=0*GradVal.V;
end
if doInputGradY
    GradVal.BY=-gradBY/SzScale;
end
if doInputGradX
    GradVal.BX=-gradBX/SzScale;
end
if ~ParamValDep.ParamStr.isFixed.Qchol
GradVal.Qchol=gradQrt;
end
%% 11.2 Parameter update

[ParamValDep,nadamStr]=BPKF_NADAM_Update(ParamValDep,nadamStr,GradVal,iBatch);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% To Do: Add Xest update
%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Population-Based Training
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if doPBT && mod(iBatch,PBTstr.SetLength)==0 && iBatch>0
    disp('Doing PBT')
    PBTcost=BPKF_PBT_Evaluate(Out,PBTstr,iBatch);
    PBTstr=BPKF_PBT_Write(PBTcost,ParamValDep,nadamStr,ModelSpec,KalSpec,PBTstr,iBatch);
    [ParamValDep,nadamStr,ModelSpec,KalSpec,PBTstr]=...
    BPKF_PBT_Exploit(PBTcost,ParamValDep,nadamStr,ModelSpec,KalSpec,PBTstr,iBatch);
    
    [Rset,Hset]=BPKF_Extract_Meas_Model(ModelSpec,numel(MeasSet));

    if strcmpi(KalSpec.BFcov(1),'s')
    Pfix=KalSpec.Pfix;
    decP=KalSpec.decP;
    decFix=KalSpec.decFix;
    elseif strcmpi(KalSpec.BFcov(1),'y')
    Pfix=KalSpec.Pfix;
    end
    if doHuber
       decHuber=ParamValDep.decHuber;
       HuberScale=ParamValDep.HuberScale;
    end
end

    
%% Adaptive KF Setup
if doAdaR
    Rbase=Rset;
    alphaR=KalSpec.alphaR;
    betaR=KalSpec.betaR;
    mixR=KalSpec.mixR;
    mixR(3)=1-(mixR(2)+mixR(1));
end
if doAdaQ
    alphaQ=KalSpec.alphaQ;
    betaQ=KalSpec.betaQ;
    mixQ=KalSpec.mixQ;
    mixQ(3)=1-(mixQ(1)+mixQ(2));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Saving Predictions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if doRecPred
if mod(iBatch,500)==1
    %% Last prediction isn't saved because not re-used in backprop
    Out.truePred{ceil(iBatch/500)}=Zrec(:,(1:(nRec-1))+ShiftBatch(1));
    if nStackEKF==1
    Out.Pred{ceil(iBatch/500)}=H*[Xrec{2:end}];
    else
        tmpRec=Uncellfun(@(xx)(xx(:,1)),Xrec(2:end));
        Out.Pred{ceil(iBatch/500)}=H*[tmpRec{:}];
    end
%    if iBatch>5000
%figure;subplot(1,2,1);plot([Xrec{:}]')
%subplot(1,2,2);plot((H*[Xrec{2:end}])','r');
%pbaspect([1 1 1])
%title(['Batch :',num2str(iBatch)])
%hold on;plot(Zrec(:,(1:(nRec-1))+ShiftBatch)','b')
%pbaspect([1 1 1])
%    title(['Batch :',num2str(iBatch)])
%    end
end
end
end
if doHuber && ~isCVrun
    usedSet=unique(abs(SetRep));
    for iH=usedSet
        usedSet0=abs(SetRep)==iH;
        med_Huber_kalSet{iH}=(1-decHuber)*HuberScale*median(med_kalHuber(:,usedSet0),2)+decHuber*med_Huber_kalSet{iH};
        med_Huber_recSet{iH}=(1-decHuber)*HuberScale*median(med_recHuber(:,usedSet0),2)+decHuber*med_Huber_recSet{iH};
    end
end
end


if doPBT
    PBTstr.FinCost=BPKF_PBT_Evaluate(Out,PBTstr,iBatch);
end

Out.FinBatch=iBatch;
Out.runTime=toc(batchTIC);
%% 12. Conclusion
[Out.Param,Out.rtQ]=BPKF_Condense_Param(ParamValDep,doInput);
if doSave
    save(ParStr.saveName,'Out')
%    delete([ParStr.saveName,'_Ahalf']);
%    delete([ParStr.saveName,'_Bhalf']);
end
Out.ParamValDep=ParamValDep;
if doAdaR
    Out.posR=posR;
    Out.Rbase=Rbase;
end
if doAdaQ
    Out.posQ=posQ;
end
if doPBT
    Out.PBT=PBTstr;
    if isfield(PBTstr,'pendingFile')&&~isempty(PBTstr.pendingFile)
        delete(PBTstr.pendingFile)
        PBTstr.pendingFile=[];
    end    
end
Out.VersStamp=VersStamp;
Out.VersNotes=VersNotes;
end