function[Out]=BPKF_Sim_DT(ooP,Start,DnSamp,tEnd,dwScale,doPlot,dt)

%% HARD-CODED for speed
W=double(ooP.Param{1});
D=double(ooP.Param{2});
C=double(ooP.Param{3});
S=double(ooP.Param{4});
V=double(ooP.Param{5});

nX=size(Start,1);
nDat=size(Start,2);

tVec=0:(dt*DnSamp):tEnd;
C=dt*C;
W=dt*W;
D=1-dt*(1-D);
nT=numel(tVec);

if size(ooP.rtQ,2)~=1
dW=diag(ooP.rtQ).*dwScale*dt;
else
dW=ooP.rtQ.*dwScale*dt;
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
    Noise=dW.*randn(nX,tEnd/dt);
    else
        %% Remark--noise isn't in same dimensions as X to improve index speed
    Noise=dW.*randn(nX,nDat,tEnd/dt);
    end
    
    
if nDat==1
    for i=1:(tEnd/dt)
%        Start=single(W*double(tanh(S.*Start+V)))+D.*Start+C+Noise(:,i);
        tt=-1+2./(1+exp(-2*(S.*Start+V)));
        Start=W*tt+D.*Start+C+Noise(:,i);
        if mod(i,DnSamp)==0
            Out(:,1+(i/DnSamp))=Start;
        end
    end
else
    for i=1:(tEnd/dt)
        tt=-1+2./(1+exp(-2*(S.*Start+V)));
        Start=W*tt+D.*Start+C+Noise(:,:,i);
        if mod(i,DnSamp)==0
            Out(:,1+(i/DnSamp),:)=Start;
        end
    end
end

else

if nDat==1
    for i=1:(tEnd/dt)
%        Start=W*tanh(S.*Start+V)+D.*Start+C;
        tt=-1+2./(1+exp(-2*(S.*Start+V)));
        Start=W*tt+D.*Start+C;
        if mod(i,DnSamp)==0
            Out(:,1+(i/DnSamp))=Start;
        end
    end
else
    for i=1:(tEnd/dt)
%        Start=W*tanh(S.*Start+V)+D.*Start+C;
        tt=-1+2./(1+exp(-2*(S.*Start+V)));
        Start=W*tt+D.*Start+C;
        if mod(i,DnSamp)==0
            Out(:,1+(i/DnSamp),:)=Start;
        end
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