function[ParamStr]=BPKF_Apply_Dependency(ParamStr,DepStr)
%% Overwrite parameter values with those used in the dependency
%% DepStr.Var.Const = vector with constant values
%% DepStr.Var.Vec [cell 1 x nGroup](nPop x nDF)
%% DepStr.Var.Scale [cell 1 x nGroup](nDF x nClass)
%% DepStr.Var.Shift [cell 1 x nGroup](1 x nClass)
%% Note--fixed values are given precedence over bounds (can violate bounds)
%% Enable the integrated representation
if nargin==1
    DepStr=ParamStr.DepStr;
    ParamStr=ParamStr.ParamStr;
end

if ~isempty(DepStr)
usedField=@(xx)(isfield(DepStr.Dep,xx)&&~isempty(DepStr.Dep.(xx)));

%% 1. Apply constant (region homogeneous) dependencies to the combined "tmpMat"
if usedField('ConstMask')||usedField('GroupMask')
nX=size(ParamStr.W,1);
%nPop=DepStr.nPop;
if usedField('ConstMask')
    tmpMat=nan(size(DepStr.Dep.ConstMask{1}));
    %% ConstMask marks the concatenated variables
for ii=1:numel(DepStr.Var.Const)
    tmpMat(DepStr.Dep.ConstMask{ii}~=0)=DepStr.Var.Const(ii);
end
else
    tmpMat=nan(size(DepStr.Dep.GroupMask{1}));
end

%% 2. Apply "group" (region heterogeneous) dependencies
if usedField('GroupMask')
    for iGroup=1:numel(DepStr.Dep.GroupMask)
        %% 2.1 Combine affine dependency model
        tmp=DepStr.Var.Shift{iGroup}+...
            DepStr.Var.Vec{iGroup}*DepStr.Var.Scale{iGroup};
        %% 2.2 Store values in tmpMat (iTy=subgroup)
        for iTy=1:size(tmp,2)
%            tmpMat(DepStr.Dep.GroupMask{iGroup}==iTy)=tmp(:,iTy);
            %% This version's faster
            tmpMat(DepStr.Dep.GroupMask{iGroup}{iTy})=tmp(:,iTy);
        end
    end
end

%figure;ScatterLine(tmpMat(1:(end/2),end-1),tmpMat(((end/2)+1):end,end-1))

%% 3. Concatenate original parameter values into matrix
tmpOrig=cell(1,numel(ParamStr.varOrder));
for ii=1:numel(ParamStr.varOrder)
tmpOrig{ii}=ParamStr.(ParamStr.varOrder{ii});
end
tmpOrig=[tmpOrig{:}];
%tmpFull=[BPKF_diag2block(tmpMat(:,1:nPop),nPop) tmpMat(:,(nPop+1):end)];

%% 4. Assign dependency-values 
tmpOrig(DepStr.isDep)=tmpMat(DepStr.isDep);

%tmpOrig=max(ParamStr.minVal,min(ParamStr.maxVal,tmpOrig));

%% 5. Store updated values
if ~strcmpi(ParamStr.varOrder{1},'W')
    error('The first parameter should always be W')
end
ParamStr.W=tmpOrig(1:nX,1:nX);
%% 5.1 If applicable, shift by base-value (depBase)
if isfield(DepStr,'depBase')
    for ii=2:numel(ParamStr.varOrder)    
    if isfield(DepStr.depBase,ParamStr.varOrder{ii})
        baseShift=DepStr.depBase.(ParamStr.varOrder{ii});
    else
        baseShift=0;
    end
    ParamStr.(ParamStr.varOrder{ii})=tmpOrig(:,nX+ii-1)+baseShift;
    end
else
    for ii=2:numel(ParamStr.varOrder)   
        ParamStr.(ParamStr.varOrder{ii})=tmpOrig(:,nX+ii-1);
    end
end

end
end

%% 6. Reformat Q from vector to matrix
if (~ParamStr.isFixed.Qchol)&&size(ParamStr.Qchol,2)==1
    ParamStr.cholQmat=diag(ParamStr.Qchol);
else
    ParamStr.cholQmat=ParamStr.cholQorig;
end


%Wtmp=ParamStr.W;
%Wmask=ParamStr.TrueMask;

%% 7 Apply parameter bounds
for ii=1:numel(ParamStr.varOrder)
    if ParamStr.doMin(ii)
       tmpName=ParamStr.varOrder{ii};
       if ParamStr.doMax(ii)
           ParamStr.(tmpName)=min(ParamStr.maxBound{ii},max(ParamStr.minBound{ii},ParamStr.(tmpName)));
       else
           ParamStr.(tmpName)=max(ParamStr.minBound{ii},ParamStr.(tmpName));
       end
    elseif ParamStr.doMax(ii)
        ParamStr.(tmpName)=min(ParamStr.maxBound{ii},ParamStr.(tmpName));
    end
end
%ParamStr.W(ParamStr.TrueMask==0)=0;
%ParamStr.W(ParamStr.TrueMask==1)=max(0,ParamStr.W(ParamStr.TrueMask==1));
%ParamStr.W(ParamStr.TrueMask==-1)=min(0,ParamStr.W(ParamStr.TrueMask==-1));

%% 7. Apply fixed (pre-specified) values where applicable
    if any(ParamStr.fixedVec,[1 2])
    for ii=find(ParamStr.fixedVec(1:numel(ParamStr.varOrder)))
        tmpName=ParamStr.varOrder{ii};
        if ParamStr.isFixed.(tmpName)
            ParamStr.(tmpName)=(ParamStr.(tmpName).*ParamStr.FreeMask{ii})+ParamStr.FixedVal{ii};
        else
            error('Fixed parameter markers are inconsistent')
        end
    end
    end


    %% Kernel overrides all other dependency components
if isfield(ParamStr,'KernParam')
    if ParamStr.ExpKern==1
for ii=1:numel(KernParam)
    ParamStr.Kern{ii}=max(0,ParamStr.Kern{ii});
    ParamStr.W(DepStr.KernMask{ii})=(ParamStr.KernSign{ii}.*ParamStr.Kern{ii}(1)).*exp(-DepStr.vecDist{ii}*max(0,ParamStr.Kern{ii}(2)));
end
    else
for ii=1:numel(KernParam)
    ParamStr.W(DepStr.KernMask{ii})=DepStr.KernFun{ii}(DepStr.vecDist{ii},ParamStr.Kern{ii});
end
    end
end