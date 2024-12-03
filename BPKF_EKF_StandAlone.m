function[Out]=BPKF_EKF_StandAlone(MeasSet,ooP,Hset,Rset,ErrMat,nStep,nStackEKF,nSim,Uinput)
%% Assumes diagonal Q
W=ooP.Param{1};
diagD=ooP.Param{2};
C=ooP.Param{3};
SS=ooP.Param{4};
V=ooP.Param{5};
rtQ=diag(ooP.rtQ);
Dmat=diag(diagD);
Q=ooP.rtQ.^2;


opts.POSDEF = true;
opts.SYM = true;

doInput=nargin>8&&~isempty(Uinput);

nX=size(W,1);

n2V=-2*V;    
n2S=-2*SS;


%probSet=@(yy)(cellfun(@(xx)(size(xx,2)-(nStep+nRec+1)),yy));
%markTrain=1:numel(MeasSet);
%probSetTrain=probSet(MeasSet(markTrain));
%probSetTrain=probSetTrain/sum(probSetTrain);

%SetRep=BPKF_Discrete_Sample(ParStr.BatchSz,probSetTrain);
%SetRep=markTrain(SetRep);


I2set=cell(size(MeasSet));
for ii=1:numel(MeasSet)
I2set{ii}=eye(size(MeasSet{ii},1));
end



%if ~isfield(ParStr,'nStack')
%    ShiftBatch=0;
%else

%% Default spacing is 5 time-steps beween members of a stack minibatch
%if isEKF && ~isempty(ParStr.nStack)
%    if isfield(ParStr,'ShiftSpace')
%        if isempty(ParStr.ShiftSpace)
%            disp('ParStr.ShiftSpace is empty, using default of 5')
%        else
%        ShiftSpace=ParStr.ShiftSpace;
%        end
%    else
%        ShiftSpace=5;
%        disp('Using default spacing of 5 within EKF stacks')
%    end
%    ShiftBatch=ShiftSpace*((1:ParStr.nStack)-1);
%else
%    ShiftBatch=0;
%end
%end
ShiftBatch=(1:nStackEKF)-1;
nStackEKF=numel(ShiftBatch);



kalSpec=ooP.Settings{2};
decFix=kalSpec.decFix;
Pfix=kalSpec.Pfix;
SimLength=kalSpec.SimLength;


simXC=cell(1,SimLength);



%% Runs initial sim twice--once to generate starting points
%% Sim from arbitrary distr. to get start points
simX=randn(nX,nSim);
  for iS=1:SimLength
    SimFun=-1+2./(1+exp(n2S.*simX+n2V));
    simX=W*SimFun+diagD.*simX+C+rtQ.*randn(nX,nSim);
  end
%% Actual sim at steady-state distr. for calculating cov  
  for iS=1:SimLength
    simXC{iS}=simX;
    SimFun=-1+2./(1+exp(n2S.*simX+n2V));
    simX=W*SimFun+diagD.*simX+C+rtQ.*randn(nX,nSim);
  end
    simXC{SimLength+1}=simX;

simMean=mean([simXC{2:end}],2);
simCov=([simXC{2:end}]-simMean)*([simXC{2:end}]-simMean)'/(SimLength*nSim-1);

Pbf=(1-decFix)*simCov+decFix*Pfix;


    %% 4.1 Measurement Specification
    BFset=cell(size(Rset));
    Pb0set=cell(size(Rset));
for iSet=1:numel(Rset)
    Hfull=Hset{iSet};
    Rtrue=Rset{iSet};
    Imat=I2set{iSet};
    Sp=linsolve(Hfull*Pbf*Hfull'+Rtrue,Imat,opts);
    BF=Pbf*Hfull'*Sp;
    Gbf=(eye(nX)-BF*Hfull);
    Pb0=Gbf*Pbf;
    BFset{iSet}=BF;
    Pb0set{iSet}=Pb0;
end

Ekal=zeros(1,nStep);
%Estore=zeros(1,nRec);

%Out=Uncellfun(@(xx)(nan),MeasSet)
Out.KalDat=cell(1,numel(MeasSet));
for iSet=1:numel(MeasSet)
MeasMark=1:nStackEKF:(size(MeasSet{iSet},2)-nStep-nStackEKF+1);
%Out00=nan(nX,nStackEKF*numel(MeasMark));
tStack=cell(1,numel(MeasMark));
Out00=cell(1,numel(MeasMark));
Emat0=ErrMat{iSet};
for iMiniBatch=1:numel(MeasMark)
    if mod(iMiniBatch,250)==1
        disp([iSet iMiniBatch numel(MeasMark)])
    end
    iB=MeasMark(iMiniBatch);
    tStack{iMiniBatch}=(1:nStackEKF)+iB+nStep-1;
    ZX0=MeasSet{iSet}(:,(iB+ShiftBatch));

    if doInput
        ZX0=ZX0-BY*Uinput{iSet}(:,(iB+ShiftBatch));
    end
    X=BFset{iSet}*ZX0;
    P=Pb0set{iSet};

    H=Hset{iSet};
    
    R=Rset{iSet};

%% 4.3 Select Measurement and Input data

Zkal=MeasSet{iSet}(:,iB+(1:(max(ShiftBatch)+(nStep))));
%Zrec=MeasSet{iSet}(:,iB+(1:(max(ShiftBatch)+(nRec)))+nStep);
if doInput
    Ukal=Uinput{iSet}(:,iB+(0:(max(ShiftBatch)+(nStep))));
%    Urec=Uinput{iSet}(:,iB+(0:(max(ShiftBatch)+(nRec)))+nStep);
    Zkal=Zkal-BY*Ukal(:,2:end);
 %   Zrec=Zrec-BY*Urec(:,2:end);
%UR=cell(1,nRec+1);
%for iRec=1:(nRec+1)
%     UR{iRec}=Urec(:,iRec+ShiftBatch);
%end
end

%PC=cell(1,nStep);



for iKal=1:nStep
 %   XC{iKal}=X;
 %   PC{iKal}=P;
 Fun=-1+(2./(1+exp((-2*SS).*X-(2*V))));
    X1t=W*Fun+diagD.*X+C;
   % if do1X
   %     Fprime=SS.*(1-Fun.^2);
   % else
    Fprime=SS.*sum(1-Fun.^2,2)/nStackEKF;
   % end
    Jac=(W.*Fprime')+Dmat;
    P1t=Jac*P*Jac'+Q;
%% 5.2 Measurement Prediction
    HP1C=H*P1t;
    yFull=Zkal(:,ShiftBatch+iKal)-Hfull*X1t;
    tmp0=Emat0*yFull;
    Ekal(iKal)=yFull(:)'*tmp0(:);

    %% 5.4 Kalman Gain
    halfS=(HP1C*H'+R)/2;
    invS=linsolve(halfS+halfS',Imat,opts);
    K=HP1C'*invS;
    %% 5.5 Correction
    X=X1t+K*yFull;
    P=P1t-K*HP1C;
    P=(P+P')/2;
end
Out00{iMiniBatch}=X;
Out.Ekal(:,iMiniBatch)=Ekal;
end
Out.Kal{iSet}=Out00;
Out.tStack{iSet}=tStack;
end
%% 6 Kalman-Free phase
%% 6.1 Predictions and Error
 %  for iRec=1:nRec
 %  Fun=-1+(2./(1+exp(-2*(SS.*X+V))));
 %  X=W*Fun+diagD.*X+C;
 %  if doInput
 %      X=X+BX*UR{iRec};
 %  end
 %  tmpE=Zrec(:,iRec+ShiftBatch)-H*X;
 %  tmp0=Emat0*tmpE;
 %  Estore(iRec)=tmpE(:)'*tmp0(:);
 %  end

end