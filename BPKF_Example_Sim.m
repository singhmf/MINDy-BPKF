
%% Small number for demo code
nRegion=20;
%% Channels-to-region ratio of .75 
%%   (or .375 in terms of pops (E+I) to chan)
nChannel=ceil(.75*nRegion);

% Number of Reps
nRep=1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% These are the main things that you might want to tinker with
%%  for real data

%% See fig 1D in the 2025 paper for explanation
    %% Number of forward-prediction steps after Kalman-filter
        nRec=5;
    %% Number of Kalman-Filtering Steps
        nStep=10;
%% Number of batches
%% This is set low for now--increase it for performance 
%% analogous the PNAS paper
NBatch=50000;
%% Number of minibatches/iteration--likewise increase this for
%% higher performance mirroring the paper
BatchSz=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Proportion of off-diagonal conns that are nonzero
%% This ~matches empirical mask
pConn=.25;

M=cell(1,nRep);
Y=cell(1,nRep);

rng shuffle


for iRep=1:nRep
disp([iRep iRep iRep iRep iRep])

%% This simulation will generate random values for:
%% brain model parameters (W,C,D,S,V)
%% Leadfield matrix (mm.H)
%% Process Noise covariance (square-root is mm.rtR)
%% Measurement Noise covariance (square-root is mm.rtQ)

%% As described in paper, randomly sampling parameters can occasionally
%% result in pathological (non-invertible) simulations, remove these here
isGood=false;maxTry=10;iTry=0;
while (~isGood)&&(iTry<maxTry)
    iTry=iTry+1;
    mm=BPKF_EI_Simulation_Example(nRegion,nChannel,pConn);
    tmp=movvar(mm.trueX,50,[],2);
    %% Criteria to drop pathalogical simulations
    isGood=mean(median(tmp,[2,3])<.01)<.5;
    if ~isGood
        disp('Failed Sim, regenerating')
    end
end
if ~isGood
   error('Too many pathological sims--something is wrong')
end

mm.Param=Uncellfun(@single,mm.Param);

%% Create mask of eligible connections
Wmask=sign(mm.Param{1});

%% Same spatial gradient for EE,EI,II,IE (up to affine transform)
%% (like in paper for data) to resolve uniqueness issues


Dep.Group=[1 -1 0 0 0 0;
           1 -1 0 0 0 0];
%% This is read as local:  Wee Wei Ce De Se Ve
%%                         Wei Wii Ci Di Si Vi
%% numbers (just 1 in this case) group factors
%% signs fix sign (i.e. Wee,Wei positive). Complex values=arbitrary sign

ModelSpec0.Dep=Dep;
%% Only interested in comparing W for these analyses
%%      with random simulations there are more concerns about weird symmetries
%%      see SI of 2025 paper for more discussion of identifiability
%% --let everything else be fixed (see paper for empirical constraints)
    ModelSpec0.fixC=single(mm.Param{3});
    ModelSpec0.fixV=single(mm.Param{5});
    ModelSpec0.fixD=single(mm.Param{2});
    ModelSpec0.fixS=single(mm.Param{4});

%% Use the mask from the simulation
ModelSpec0.Wmask=Wmask;
%% using tanh nonlinearity
ModelSpec0.SigTy='tanh';
%% Not fitting Q this time
ModelSpec0.freeQ='n';
%% Using simulation process-noise covariance
ModelSpec0.Q=single(mm.rtQ*mm.rtQ');
%% Measurement model has leadfield (mm.H) for excitatory, zeros for inhibitory
ModelSpec0.H={single([mm.H zeros(nChannel,nRegion)])};
%% Using simulation measurement noise covariance
ModelSpec0.R={single(mm.rtR*mm.rtR')};

%% Using EKF
kalSpec.KFtype='E';
%% Initial-covariance estimate done by simulation 
kalSpec.BFcov='s';
%% Specifying some simulation factors
kalSpec.SimLength=40;
kalSpec.nSim=20;
kalSpec.nSaveStart=250;

%% Specifying the length of Kalman-Filtering
kalSpec.nStep=nStep;
%% Length of free-running predictions (no KF)
kalSpec.nRec=nRec;

%% AR-smoothing of initial covariance estimates
%%      more efficient than running long sims each iteration
kalSpec.Pfix=single(eye(2*nRegion)/4);
kalSpec.decP=.95;
kalSpec.decFix=.2;
kalSpec.Pbase=single(eye(2*nRegion)/4);

%% Start calculating errors starting at start of KF (default)
kalSpec.minStep=0;

%% Number of minibatches/iteration
ParStr0.BatchSz=BatchSz;
%% Number of total iterations
ParStr0.NBatch=NBatch;
%% No penalty for steady-state error estimation
ParStr0.CovErr=0;
%% Optionally Record W-estimates every 100 iterations
    %ParStr0.recW='y';
    %ParStr0.recSpace=100;
%% 80 time-samples per minibatch
ParStr0.nStack=80;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gradient Parameterization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gradient optimizer hyperparameters (NADAM by default)
GradSpec.reg=.0001;
GradSpec.dec1=.98;
GradSpec.dec2=.99;
GradSpec.Rate=.0001;
%% Gradient clipping: relative magnitude
GradSpec.Clip=3;
%% Moving-average estimate of the magnitude to clip gradients (last k-steps)
GradSpec.ClipLength=200;
GradSpec.ClipType='Moving';


%%%%%%%%%%%%%%%%%%
%% Formatting Data
%%%%%%%%%%%%%%%%%%

Xguess=cell(1,size(mm.trueX,3));
MeasSet=cell(1,size(mm.trueX,3));
for ii=1:size(mm.trueX,3)
   %% Xguess is deprecated--fill with garbage0
        Xguess{ii}=single(randn(size(mm.trueX,[1 2]))/100);
   %% Each element of meas-set is a simulated recording 
        MeasSet{ii}=single(mm.meas(:,:,ii));
end

%% Fit model
Y{iRep}=BPKF_Full(Xguess,MeasSet,ParStr0,kalSpec,ModelSpec0,GradSpec);

Wbase=(Wmask(1:nRegion,1:nRegion).*(1-eye(nRegion)));

predWEE=Y{iRep}.Param{1}(1:nRegion,1:nRegion);
predWEI=Y{iRep}.Param{1}((1+nRegion):(2*nRegion),1:nRegion);


trueWEE=mm.Param{1}(1:nRegion,1:nRegion);
trueWEI=mm.Param{1}((1+nRegion):(2*nRegion),1:nRegion);

figure;
subplot(2,2,1);
ScatterLine(predWEE(Wbase~=0),trueWEE(Wbase~=0));
pbaspect([1 1 1])
subplot(2,2,2);
ScatterLine(predWEI(Wbase~=0),trueWEI(Wbase~=0));
pbaspect([1 1 1])
end
