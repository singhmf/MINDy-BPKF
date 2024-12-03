function[Ymeas,OutStruct]=BPKF_PreProc_Proj(Ymeas,Hfull,EigThresh,R)
%% Rank-reducing projection
%% Ymeas=data
%% Hfull=Parcellated LF
%% EigThresh= Relative threshold for censoring SVs of H (0<=EigThresh<1)
%% R=nChan x nChan noise covariance

if EigThresh>=1
    error('EigThresh should be <1')
end
if ~iscell(Ymeas)
    [Hforw,Hmeas]=Tmp_EigFun(Hfull,EigThresh);
    origCov=cov(Ymeas');
    dropDim=eye(size(Hfull,1))-Hmeas'*Hmeas;
    dropCov=dropDim*origCov*dropDim';
    Ymeas=Hmeas*Ymeas;
  %  Hnew=Hforw;
    %% Assume that the data is normal iid--rescale based upon std
    Yscale=std(Ymeas(:));%/sqrt(mean(sum(Hforw.^2,2),1));
    Ymeas=Ymeas/Yscale;
    Hnew=Hforw/sqrt(mean(sum(Hforw.^2,2),1));
    Rnew=Hmeas*R*Hmeas'/(Yscale^2);
    dropCov=dropCov/(Yscale^2);
    ProjMeas=Hmeas;
else
    Yscale=zeros(size(Ymeas));
    dropCov=cell(size(Ymeas));
    origCov=cell(size(Ymeas));
    if ~iscell(Hfull)
        Hfull=repmat({Hfull},size(Ymeas));
        wasHcell=false;
    else
        ProjMeas=cell(size(Hfull));
        Hnew=cell(size(Hfull));
        wasHcell=true;
    end
    if ~iscell(R)
        R=repmat({R},size(Ymeas));
        wasRcell=false;
    else
        Rnew=cell(size(R));
        wasRcell=true;
    end
    for ii=1:numel(Ymeas)
    [Hforw,Hmeas]=Tmp_EigFun(Hfull{ii},EigThresh);
    origCov{ii}=cov(Ymeas{ii}');    
    dropDim=eye(size(Hfull{ii},1))-Hmeas'*Hmeas;
    Ymeas{ii}=Hmeas*Ymeas{ii};
        %% Assume that the data is normal iid--rescale based upon std
    Yscale(ii)=std(Ymeas{ii}(:));
    Hforw=Hforw/sqrt(mean(sum(Hforw.^2,2),1));%/sqrt(mean(sum(Hforw.^2,2),1));
    dropCov{ii}=dropDim*origCov{ii}*dropDim'/(Yscale(ii)^2);
    Ymeas{ii}=Ymeas{ii}/Yscale(ii);
    if wasRcell
    Rnew{ii}=Hmeas*R{ii}*Hmeas';
    elseif and(~wasRcell,wasHcell)
        error('Number of H and R cases must agree')
    else
        Rnew=Hmeas*R{ii}*Hmeas'/(Yscale(ii)^2);
    end
    if wasHcell
        Hnew{ii}=Hforw;
        ProjMeas{ii}=Hmeas;
    else
        Hnew=Hforw;
        ProjMeas=Hmeas;
    end

    end
end
OutStruct.Hnew=Hnew;
OutStruct.ProjMeas=ProjMeas;
OutStruct.Rnew=Rnew;
OutStruct.Yscale=Yscale;
OutStruct.dropCov=dropCov;
OutStruct.origCov=origCov;
OutStruct.note='Yscale has already been factored into Hnew, Rnew, and dropCov';
end







function[Hx,Hy]=Tmp_EigFun(H0,E0)
%% Hx: remaining H matrix
%% Hy: what to multiply Ymeas by
[u,s,v]=svd(H0); 
Spass=(sum(s,2)/max(s(:)))>=E0;
Hy=u(:,Spass)';
Hx=s(Spass,:)*v';
end

