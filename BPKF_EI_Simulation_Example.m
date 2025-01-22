function[ooP]=BPKF_EI_Simulation_Example(nX,nH,pConn,rtQ)
%%  Remark: several parameters in simulations are fixed (e.g. S)
%%  This is b/c (via a change-of-variables) their distribution would be absorbed
%%      in another variable (which isn't fixed).
%%      Better to reduce to the minimum unique dimensions re:dynamical-conjugacy
%%      for hyper-distribution setting instead of lots of redundant distributions

%%      i.e. x(t+1)=W0 psi(s0*x+v0)+D+c0+rtQ0*epsilon
%%           y(t+1)=H0x+rtR*eta

%%  is isomorphic to   q<-s0x+v0    z<-inv(rtR0)*(y+H0*v0)
%%           q(t+1)=Wpsi(q)+D+c+rtQ*epsilon
%%           z(t+1)=Hq+eta
%%           W=diag(s0)W0,  c=s0*c0+(I-D)v0,   rtQ=s0*rtQ0 H=inv(rtR0)*H0*inv(s0)
%%           Hence, no need to specify distributions for R, V, S in simulations


nRep=20;
if nargin<3 || isempty(pConn)
    disp('Using default pConn')
    pConn=25;
end
if pConn<1 && pConn~=0
    disp('Assuming pConn entered as decimal--converting to percent')
    pConn=100*pConn;
end

tRun=2000;
dropStart=1;

%% Random, symmetric base mask
W0=rand(nX);W0=W0+W0';
tmp=prctile(OffDiag(W0),100-pConn);
Wmask=eye(nX)+(NoDiag(W0)>=tmp);

%% Generate connection values independent of mask, from hyper-distribution
%% (hence asymmetric, some will still be small)
    %% Random low-rank structure
    nRank=ceil(nX/4);
    tmp1=(rand(nX,nRank).^3+rand(nX,nRank)/5);
    tmp2=(rand(nX,nRank).^3+rand(nX,nRank)/5);
    WE=NoDiag(Wmask).*(tmp1*tmp2')/max(tmp1*tmp2',[],'all')*4/nX;
    WE=4*WE/3;
    %% Add to random sparse structure
    WE=NoDiag(Wmask).*(WE+(16/nX)*rand(nX).^3);

%% Same thing but for the EI connections
    tmp1=(rand(nX,nRank).^3+rand(nX,nRank)/5);
    tmp2=(rand(nX,nRank).^3+rand(nX,nRank)/5);
    WI=NoDiag(Wmask).*(tmp1*tmp2');WI=(WI/max(WI,[],'all'))*4/nX;
    WI=NoDiag(Wmask).*(WI+(16/nX)*rand(nX).^3);

%% Base ratio of local strengths (multiplied by random vector in next few lines)
RecBase=[2 -5;5 -1]/5;
%% randomly generate spatial pattern
recVecE=abs(1+.1*rand(nX,1).^3);
%% same "pattern" for E and I
recVecI=recVecE;

%% Rescaling of random pattern
recE=[diag(recVecE.*RecBase(1,1));diag(recVecE.*RecBase(2,1))];
recI=[diag(recVecI.*RecBase(1,2));diag(recVecI.*RecBase(2,2))];

Wfull=[recE recI];
Wfull=Wfull+([WE zeros(nX);WI zeros(nX)]);

%% Fixed shift-terms at zero, to reduce number of pathological (saturated) cases
%% Users can change to a small random vector if desired...
C=0;

%% V is absorbed in C distribution via a change-of-variables...
V=0;

%% Random scalar values for decay (since scalar D used in empirical model)
DE=max(.5,min(.9,.65+rand(1)/50));
DI=max(.5,min(.9,.8+rand(1)/50));

%% Fixed, but no loss-of-generality as random S is absorbed in the W distribution
%% via a change-of-variables
SE=2.5;
SI=2.5;

%% Fully random lead-field
H=randn(nH,nX);
%% Random uncorrelated measurement noise but no loss of generality
    %% (correlated cases behave isomorphic to rotating H by PCs of R)
rtR=diag(.2+rand(nH,1)/10)+rand(1)/10;

%% Scalar values to vector
D=repelem([DE;DI],nX,1);
S=repelem([SE;SI],nX,1);

ooP.Param={Wfull,D,C,S,V};
ooP.H=H;

if nargin<4 || isempty(rtQ)
        disp('Using default rtQ')
    ooP.rtQ=diag(repelem([.5;.5],nX,1));
elseif numel(rtQ)==1 %#ok<ISCL>
    ooP.rtQ=diag(repelem(rtQ,2*nX,1));
elseif numel(rtQ)==2
    ooP.rtQ=diag(repelem(rtQ(:),nX,1));
end
ooP.rtR=rtR;

%% Start simulation to get up to long-term distribution
tmpX=BPKF_Sim(ooP,randn(2*nX,nRep),1,tRun,ooP.rtQ(1),'n');
tmpX=tmpX(:,(1+dropStart):end);

%% Data for training drawn from long-term distribution (more realistic)
ooP.trueX=BPKF_Sim(ooP,randn(2*nX,nRep).*std(tmpX(:,:),[],2),1,tRun,1,'n');


%% Make measurements by adding noise and multiplying by leadfield (H)
ooP.meas=zeros(nH,size(ooP.trueX,2),nRep);
ooP.meas(:,:)=rtR*randn(nH,nRep*size(ooP.trueX,2))+H*ooP.trueX(1:nX,:);
end

