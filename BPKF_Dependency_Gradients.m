function[GradVal]=BPKF_Dependency_Gradients(ParamStr,DepStr,GradVal)
%% Overwrite parameter values with those used in the dependency
%% DepStr.Var.Const = vector with constant values
%% DepStr.Var.Vec [cell 1 x nGroup](nPop x nDF)
%% DepStr.Var.Scale [cell 1 x nGroup](nDF x nClass)
%% DepStr.Var.Shift [cell 1 x nGroup](1 x nClass)

%% NOTE: This version just generates the gradients for dependent variables--
%%      leaves the existing (independent) gradients (to be overwritten)


%% Enable the integrated representation (alternative input: ParamValDep+GradVal)
if nargin==2
    GradVal=DepStr;
    DepStr=ParamStr.DepStr;
    ParamStr=ParamStr.ParamStr;
end

usedField=@(xx)(isfield(DepStr.Dep,xx)&&~isempty(DepStr.Dep.(xx)));

gradOrig=cell(1,numel(ParamStr.varOrder));
for ii=1:numel(ParamStr.varOrder)
gradOrig{ii}=GradVal.(ParamStr.varOrder{ii});
end
gradOrig=[gradOrig{:}];

%if usedField('ConstMask')
%if usedField('ConstMask')||usedField('GroupMask')
%nX=size(ParamStr.W,1);
%nPop=DepStr.nPop;

if usedField('ConstMask')
    for ii=1:numel(DepStr.Var.Const)
        GradVal.Const(1,ii)=sum(gradOrig(DepStr.Dep.ConstMask{ii}~=0),[1 2]);
    end
end
if usedField('GroupMask')
    for iGroup=1:numel(DepStr.Dep.GroupMask)
    tmpMat=zeros(sum(DepStr.Dep.GroupMask{iGroup}{1},[1 2]),numel(DepStr.Dep.GroupMask{iGroup}));
    for iClass=1:numel(DepStr.Dep.GroupMask{iGroup})
        tmpMat(:,iClass)=gradOrig(DepStr.Dep.GroupMask{iGroup}{iClass});
    end
    GradVal.Shift{iGroup}=sum(tmpMat,1);
    GradVal.Vec{iGroup}=tmpMat*DepStr.Var.Scale{iGroup}';
    GradVal.Scale{iGroup}=DepStr.Var.Vec{iGroup}'*tmpMat;
    end
end
if usedField('KernDist')
    if ParamStr.ExpKern==1
        for ii=1:numel(DepStr.KernDist)
        dKern=ParamStr.KernSign{ii}.*GradVal.W(DepStr.KernMask{ii}).*exp(-DepStr.vecDist*ParamStr.KernParam{ii}(2));
        GradVal.Kern{1,ii}(1)=sum(dKern,'all');
        GradVal.Kern{1,ii}(2)=-ParamStr.KernParam{ii}(1)*sum(dKern.*DepStr.vecDist,'all');
        end
    else
        for ii=1:numel(DepStr.KernDist)
            dKern=GradVal.W(DepStr.KernMask{ii});
            GradVal.Kern{1,ii}=DepStr.gradFun{ii}(dKern,DepStr.vecDist,ParamStr.KernParam{ii});
        end
    end
end
end

