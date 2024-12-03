function[DepStr,ParamStr,nadamStr]=BPKF_Initialize_Dependency(ModelSpec,ParamStr,nadamStr)
%% For now, only contain functionality for diagonal Q!!

%% Constant dependencies=same value for all elements within group:
%%   block array of indices by subpop with sign or complex (for free sign)
%%   0=no constant dependency

%% Unique fct but reporting at most 1 nan
unique1nan=@(xx)([unique(xx(~isnan(xx)));xx(find(isnan(xx),1,'first'))]);

nX=size(ParamStr.W,1);

%% 1 Initialize Settings
if isfield(ModelSpec,'Dep')&&(~isempty(ModelSpec.Dep))
DepStr=ModelSpec.Dep;
RandScale=ModelSpec.RandScale;
usedDep=@(xx)(isfield(DepStr,xx)&&~isempty(DepStr.(xx)));

%% 1.1 Set up degrees of freedom for each dependency
if usedDep('Group')&&~usedDep('dfGroup')
    disp('Assuming 1 df per group')
    DepStr.dfGroup=ones(1,numel(unique(DepStr.Group(DepStr.Group~=0))));
elseif usedDep('Group')&& and(numel(DepStr.dfGroup)==1,numel(unique(DepStr.Group(DepStr.Group~=0)))>1)
    disp('Assuming same df per group')
    DepStr.dfGroup=repmat(DepStr.dfGroup,1,numel(unique(DepStr.Group(DepStr.Group~=0))));
end

%% 2. Ensure consistency
if usedDep('Const')
if usedDep('Group')
   szGroup=size(DepStr.Group);
   szConst=size(DepStr.Const);
   tmpConst=DepStr.Const;
   tmpGroup=DepStr.Group;
   %% 2.1 If using implied dimensions for group or constant, resize accordingly 
   if szGroup(1)>szConst(1)
       tmpConst=[repmat(tmpConst(1),szGroup(1),szGroup(1)) repmat(tmpConst(2:end),szGroup(1),1)];
       disp('Resizing Const to match Group notation')
   elseif szGroup(1)<szConst(1)
       tmpGroup=[repmat(tmpGroup(1),szConst(1),szConst(1)) repmat(tmpGroup(2:end),szConst(1),1)];
       disp('Resizing Group to match Const notation')
   end
   %% 2.2 Check for inconsistencies/overlap
   if any(size(tmpConst)~=size(tmpGroup))
       error(['Incompatable array sizes for Const (',num2str(size(DepStr.Const)),...
           ') and Group (',num2str(size(DepStr.Group)),') dependencies']);
   end

   %% 2.3 Determine which variables have dependencies (isDep)
   isDep=(tmpConst~=0)+(tmpGroup~=0);
   if any(isDep>1,[1 2])
       error('Overlapping Const and Group dependencies')
   end
else
    isDep=DepStr.Const~=0;
end
elseif usedDep('Group')
    isDep=DepStr.Group~=0;
else
    isDep=false;
end

if numel(isDep)==1
    DepStr.isDep=isDep;
    DepStr.nPop=1;
else
    %% Expand dependencies to full combined-parameter matrix
    nPop=size(isDep,1);
    DepStr.isDep=[BPKF_diag2block(repelem(isDep(:,1:nPop),nX/nPop,1),nPop)==1 repelem(isDep(:,(nPop+1):end)==1,nX/nPop,1)];
end


%% 3 Setup constant dependencies

if usedDep('Const')
    if any(isnan(DepStr.Const),[1 2])
        error('Use i not nan to indicate sign-free (enables group marking 1i, 2i, etc.)')
    end
    nPop=size(DepStr.Const,1);
    %% This error can be removed later
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% REMARK: Can't include BX, BY matrices in dependency
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if size(DepStr.Const,2)~=(nPop+numel(ParamStr.varOrder)-1)
        disp(strcat('Dependency constants specification is mis-sized'))
            error(strcat('Dependency constants should be [',ParamStr.varOrder,']'));
    end
    %% 3.1 Switch to block format
    ConstMask=[kron(DepStr.Const(:,1:nPop),eye(nX/nPop)),...
        repelem(DepStr.Const(:,(nPop+1):end),nX/nPop,1)];
    %% Signmat stores as 0,+/-1,nan instead of using complex numbers
    ConstSign=sign(real(DepStr.Const));
    ConstSign(imag(DepStr.Const)~=0)=nan;
   
    ConstGroups=unique(abs(DepStr.Const(DepStr.Const~=0)))';
    %% Make sure numbers aren't skipped
    if numel(ConstGroups)~=max(abs(DepStr.Const(:)))
        error(['Skipped an index in Const dependencies: ',num2str(ConstGroups)])
    end
    %% This only occurs when there're e.g. 1 & -1, or 1 & 1i, etc
    %% constant dependencies are forced to be exactly the same between different variables (can't have different signs)
    if numel(ConstGroups)~=numel(unique(DepStr.Const(DepStr.Const~=0)))
           error('Mixed sign assignments for elements with the same constant-dependency')
    end
    DepStr.Dep.ConstMarkers=ConstGroups;
    %% 3.2 Setup NADAM initialization
    nadamStr.mConst=zeros(1,numel(ConstGroups));
    nadamStr.nConst=zeros(1,numel(ConstGroups));
    
    for ii=1:numel(ConstGroups)
       DepStr.Dep.ConstMask{ii}=abs(ConstMask)==ConstGroups(ii);
       depSign=unique1nan(ConstSign(abs(DepStr.Const)==ii));
       if numel(depSign)~=1
           error('Mixed sign assignments for elements with the same constant-dependency')
       else
           DepStr.Dep.ConstSign(1,ii)=depSign;
       end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %% 3.3 Set Initial random values for constants
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       if isnan(depSign)
           DepStr.Var.Const(ii)=RandScale*randn(1);
       else
           DepStr.Var.Const(ii)=depSign*abs(randn(1))*RandScale;
       end
    end
end


if usedDep('Group')
    if any(isnan(DepStr.Group),[1 2])
        error('Use i not nan to indicate sign-free (enables group marking 1i, 2i, etc.)')
    end
    nPop=size(DepStr.Group,1);
    if size(DepStr.Group,2)~=(nPop+numel(ParamStr.varOrder)-1)
        ErrTmp=[ParamStr.varOrder;repmat({','},size(ParamStr.varOrder))];
        error(strcat('Dependency groups specification is mis-sized: ',...
            ' should be [',[ErrTmp{:}],']'));
    end
    %% Signmat stores as 0,+/-1,nan instead of using complex numbers
    GroupSign=sign(real(DepStr.Group));
    GroupSign(imag(DepStr.Group)~=0)=nan;
   % Out.Dep.GroupSign=GroupSign;
    GroupIndx=unique(abs(DepStr.Group(DepStr.Group~=0)));
    GroupIndx=reshape(GroupIndx,1,numel(GroupIndx));
        %% Make sure numbers aren't skipped
    if numel(GroupIndx)~=max(abs(DepStr.Group(:)))
        error(['Skipped an index in Group dependencies: ',num2str(GroupIndx)])
    end
    DepStr.Dep.GroupMarkers=GroupIndx;
%    nadamStr.mShift=zeros(1,numel(GroupIndx));
%    nadamStr.nShift=zeros(1,numel(GroupIndx));

    for iGroup=GroupIndx
        groupMatch=find(abs(DepStr.Group)==iGroup);
        tmpMask=zeros(size(DepStr.Group));
        for iClass=1:numel(groupMatch)
            tmpMask(groupMatch(iClass))=iClass;
        end        
        nadamStr.mScale{1,iGroup}=zeros(DepStr.dfGroup(iGroup),numel(groupMatch));
        nadamStr.nScale{1,iGroup}=zeros(DepStr.dfGroup(iGroup),numel(groupMatch));
        nadamStr.mShift{1,iGroup}=zeros(1,numel(groupMatch));
        nadamStr.nShift{1,iGroup}=zeros(1,numel(groupMatch));
        
        nadamStr.mVec{1,iGroup}=zeros(nX/nPop,DepStr.dfGroup(iGroup));
        nadamStr.nVec{1,iGroup}=zeros(nX/nPop,DepStr.dfGroup(iGroup));
        %% Can use DepVecGroups as a marker for converting concatenated arrays back to cells
        nadamStr.Groups.Scale{1,iGroup}=iGroup*ones(1,numel(groupMatch));
        nadamStr.Groups.Shift{1,iGroup}=iGroup*ones(1,numel(groupMatch));
        nadamStr.Groups.Vec{1,iGroup}=iGroup*ones(1,DepStr.dfGroup(iGroup));
            %% switch to block format
        GroupMask=[kron(tmpMask(:,1:nPop),eye(nX/nPop)),...
            repelem(tmpMask(:,(nPop+1):end),nX/nPop,1)];
     
        
        %Out.Dep.GroupMask{1,iGroup}=GroupMask;
        %% This version's faster to apply
        for iClass=1:numel(groupMatch)
            DepStr.Dep.GroupMask{1,iGroup}{iClass}=GroupMask==iClass;
        end


        
        tmpVec=abs(randn(nX/nPop,DepStr.dfGroup(iGroup)));
        DepStr.Var.Vec{1,iGroup}=tmpVec./max(tmpVec,[],1);
        tmpScale=RandScale*randn(DepStr.dfGroup(iGroup),numel(groupMatch));
        tmpShift=RandScale*randn(1,numel(groupMatch));
        for iClass=1:numel(groupMatch)
            DepStr.Dep.ScaleSign{iGroup}(1,iClass)=GroupSign(groupMatch(iClass));
            DepStr.Dep.ShiftSign{iGroup}(1,iClass)=GroupSign(groupMatch(iClass));
            if isnan(GroupSign(groupMatch(iClass)))
                DepStr.Var.Scale{1,iGroup}(:,iClass)=tmpScale(:,iClass);
                DepStr.Var.Shift{1,iGroup}(1,iClass)=tmpShift(iClass);
            else
                DepStr.Var.Scale{1,iGroup}(:,iClass)=GroupSign(groupMatch(iClass))*abs(tmpScale(:,iClass));
                DepStr.Var.Shift{1,iGroup}(1,iClass)=GroupSign(groupMatch(iClass))*abs(tmpShift(iClass));
            end
        end
    end
end
if usedDep(KernDist)
    if ~iscell(KernDist)
        DepStr.KernDist={DepStr.KernDist};
    end
    nKern=numel(DepStr.KernDist);
    if isfield(DepStr,'ExpKern')&&DepStr.ExpKern==1
        DepStr.dimKernParam=2;
    end
    if ~iscell(DepStr.KernFun)
        DepStr.KernFun={DepStr.KernFun};
        DepStr.gradKern={DepStr.gradKern};
    end
    %% KernFun should take arguments(Distance, Parameters) and return in the same dim. as Distance
    %% gradKern should take arguments(dW, Distance, Parameters) and return in the dimension Parameter
    if numel(DepStr.KernFun)==1
        %% Function that defines Kernel
        DepStr.KernFun=repelem(DepStr.KernFun,1,nKern);
        %% Analytic gradient of Kernel w.r.t. parameters
        DepStr.gradKern=repelem(DepStr.gradKern,1,nKern);
    end
    %% vectorized distance for unmasked entries (so don't need to recalculate)
    for ii=1:numel(DepStr.KernDist)
         %% KernMask is signed
         DepStr.KernSign{ii}=DepStr.KernMask{ii}(DepStr.KernMask~=0);
         %% Condense the sign into a single element if homogeneous
         if numel(unique(DepStr.KernSign{ii}))==1
             DepStr.KernSign{ii}=DepStr.KernSign{ii}(1);
         end
         DepStr.KernMask{ii}=DepStr.KernMask{ii}~=0;
         DepStr.vecDist{1,ii}=DepStr.KernDist{ii}(DepStr.KernMask{ii});
         DepStr.vecDist{1,ii}=DepStr.vecDist{1,ii}(:);
         szKern=DepStr.dimKernParam;
         ParamStr.KernParam{1,ii}=abs(rand(szKern));
         nadamStr.mKern{1,ii}=zeros(szKern);
         nadamStr.nKern{1,ii}=zeros(szKern);
    end
end
else
    DepStr=[];
    disp('No dependencies entered')
end
end