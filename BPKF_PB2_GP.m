function[Out,predE,predSD]=BPKF_PB2_GP(X,Y,tInd,Xpend,tPend,Ep,alph,lB,uB,Tnow,GPR_tLength,doMultiLength)

tInd=reshape(tInd,numel(tInd),1);
if isempty(tPend)
    tPend=[];
    tPend=reshape(tPend,1,numel(tPend));
end

if isempty(Xpend)
    Xpend=[];
end

nSamp=size(X,1);
nX=size(X,2);
nFull=nSamp+size(Xpend,1);


%% Center
if numel(Y)>3
yMean=mean(Y);
else
    yMean=0;
end
Y=Y-yMean;

%% Test points to initialize fmincon
nTest=25;
%% Scaling ratio to sample unbounded variables
pRatio=1.5;
nPt=5;


%% Only use nearby time-points in estimating variances
t1=abs(tInd-Tnow)<=GPR_tLength;

if doMultiLength
    KernFct='ardsquaredexponential';
else
    KernFct='squaredexponential';
end
if isempty(GPR_tLength)
    Mdl=fitrgp(X,Y,'KernelFunction',KernFct,'BasisFunction','none');
else
    Mdl=fitrgp(X(t1,:),Y(t1),'KernelFunction',KernFct,'BasisFunction','none');
end

tmp=Mdl.KernelInformation.KernelParameters;

sigmaL=tmp(1:(end-1))';
sigmaF=tmp(end);
sigmaNoise=Mdl.Sigma^2;



Xfull=[X;Xpend];

rtL=sqrt(2)*sigmaL;
w_Xfull=Xfull./rtL;
w_dX=X./rtL;

distFull=sum(((permute(w_Xfull,[1 3 2])-permute(w_Xfull,[3 1 2])).^2),3);
Xdist=distFull(1:nSamp,1:nSamp);

tWeight_Pend=sqrt(1-Ep).^abs([tInd;tPend]-[tInd;tPend]');
tWeight=tWeight_Pend(1:nSamp,1:nSamp);

Kmat=sigmaF^2*exp(-Xdist).*tWeight;
KmatFull=sigmaF^2*exp(-distFull).*tWeight_Pend;

tK=sqrt(1-Ep).^abs(Tnow-tInd);
tK_P=sqrt(1-Ep).^abs((Tnow+1-[tInd;tPend]));

ktFun=@(xx)(sigmaF^2*tK.*exp(-sum((xx./rtL-w_dX).^2,2)));

ktFunP=@(xx)(sigmaF^2*tK_P.*exp(-sum((xx./rtL-w_Xfull).^2,2)));

invCovY=(Kmat+sigmaNoise^2*eye(nSamp))\Y;
invCovP=inv(KmatFull+sigmaNoise^2*eye(nFull));

sigmaF2=sigmaF^2;

if alph==0
fCost=@(xx)(ktFun(xx)'*invCovY);
else
fCost=@(xx)(ktFun(xx)'*invCovY-alph*sqrt(sigmaF2-ktFunP(xx)'*invCovP*ktFunP(xx))); %#ok<MINV> 
end

oo=optimoptions('fmincon');
oo.Display='off';


%% Test points to determine initialization

xMax=max(X,[],1);
xMin=min(X,[],1);
cMean=mean(X,1);
lb0=cMean-pRatio*(cMean-xMin);
ub0=cMean+pRatio*(xMax-cMean);

if isempty(lB)
    lB=lb0;
end
if isempty(uB)
    uB=ub0;
end
lB(isinf(lB))=lb0(isinf(lB));
uB(isinf(uB))=ub0(isinf(uB));


xTest=lB+rand(nTest,nX).*(uB-lB);
startCost=zeros(nTest,1);

for ii=1:nTest
    startCost(ii)=fCost(xTest(ii,:));
end
[~,xM]=mink(startCost,nPt);
x0=xTest(xM,:);

Out=zeros(nPt,nX);
cOut=zeros(nPt,1);
for iT=1:nPt
[Out(iT,:),cOut(iT)]=fmincon(fCost,x0(iT,:),[],[],[],[],lB,uB,[],oo);
end
[~,minMark]=min(cOut);
%disp([Out cOut])
Out=Out(minMark,:);
%[Out2,xe2]=fmincon(fCost,xTest(1,:),[],[],[],[],lb,ub,[],oo);


predE=yMean+ktFun(Out)'*invCovY;
predSD=sqrt(sigmaF2-ktFunP(Out)'*invCovP*ktFunP(Out)); %#ok<MINV> 

end






