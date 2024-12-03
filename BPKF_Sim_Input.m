function[Out]=BPKF_Sim_Input(ooP,Start,DnSamp,tEnd,dwScale,funcInput,doPlot)
%% funcInput can be an array or a function handle that accepts time vectors


%% HARD-CODED for speed
W=ooP.Param{1};
D=ooP.Param{2};
C=ooP.Param{3};
S=ooP.Param{4};
V=ooP.Param{5};

nX=size(Start,1);
nDat=size(Start,2);
dt=1;
tVec=0:(dt*DnSamp):tEnd;
C=dt*C;
W=dt*W;
D=1-dt*(1-D);
nT=numel(tVec);

dW=diag(ooP.rtQ)*dwScale*dt;

if isa(funcInput,'function_handle')
    funcInput=funcInput(dt*(1:(tEnd/dt)));
end


if nDat==1
    Out=nan(nX,nT);
    Out(:,1)=Start;
else
    Out=nan(nX,nT,nDat);
    Out(:,1,:)=Start;
end
if mean(dW(:)~=0)~=0
    if nDat==1
    NetInput=dW.*randn(nX,tEnd/dt)+funcInput;
    else
    %% Remark--noise isn't in same dimensions as X to improve index speed
    NetInput=dW.*randn(nX,nDat,tEnd/dt)+permute(funcInput,[1 3 2]);
    end
else
    if nDat==1
    NetInput=funcInput;
    else
    %% Remark--noise isn't in same dimensions as X to improve index speed
    NetInput=permute(funcInput,[1 3 2]);
    end
end
    
    
if nDat==1
    for i=1:(tEnd/dt)
        Start=single(W*double(tanh(S.*Start+V)))+D.*Start+C+NetInput(:,i);
        if mod(i,DnSamp)==0
            Out(:,1+(i/DnSamp))=Start;
        end
    end
else
    for i=1:(tEnd/dt)
        Start=W*tanh(S.*Start+V)+D.*Start+C+NetInput(:,:,i);
        if mod(i,DnSamp)==0
            Out(:,1+(i/DnSamp),:)=Start;
        end
    end
end


if ~isempty(doPlot)&&strcmpi(doPlot(1),'y')
    figure
    if nDat==1
        plot(tVec,Out);
    else
        plot(tVec,squeeze(Out(:,:,1)));
    end
end
end