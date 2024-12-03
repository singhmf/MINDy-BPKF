function[ParamValDep,nadamStr]=BPKF_NADAM_Update(ParamValDep,nadamStr,GradVal,iBatch)
%% Overwrite parameter values with those used in the dependency
%% DepStr.Var.Const = vector with constant values
%% DepStr.Var.Vec [cell 1 x nGroup](nPop x nDF)
%% DepStr.Var.Scale [cell 1 x nGroup](nDF x nClass)
%% DepStr.Var.Shift [cell 1 x nGroup](1 x nClass)

%% Remark: This version doesn't allow sign-reversals for patterns that are signed (i.e., 1-vec)

%% FOR SIMPLICITY: CLIP BEFORE APPLYING DEPENDENCY

clipInnovation=.2;

vNames=ParamValDep.ParamStr.varOrder;

if ParamValDep.ParamStr.doReg
    %% Make sure parameters already have dependencies evaluated
    ParamValDep.ParamStr=BPKF_Apply_Dependency(ParamValDep);
    if ParamValDep.Reg.doMINDy
        %% WKmask further 
        WKgrad=GradVal.W.*ParamValDep.Reg.MINDy.WKmask;
    end
end

if ~isempty(ParamValDep.ParamStr.dropGrad)
for ii=1:numel(vNames)
    if isfield(ParamValDep.ParamStr.dropGrad,vNames{ii})
        GradVal.(ParamValDep.ParamStr.varOrder{ii})=GradVal.(ParamValDep.ParamStr.varOrder{ii}).*...
            ParamValDep.ParamStr.dropGrad.(vNames{ii});
    end
end
end




if ParamValDep.ParamStr.doReg
    %% Make sure parameters already have dependencies evaluated
    ParamValDep.ParamStr=BPKF_Apply_Dependency(ParamValDep);
    if ParamValDep.Reg.doMINDy
        Mdat=ParamValDep.Reg.MINDy;
        Wk1=ParamValDep.ParamStr.Wk1;
        Wk2=ParamValDep.ParamStr.Wk2;
        %% Assuming single MINDy block for all of W
       if isempty(Mdat.PreBlock)
          GradVal.Wk1{1}=WKgrad*Wk2{1}'-Mdat.Sp1.*sign(Wk1{1});
          GradVal.Wk2{1}=Wk1{1}'*WKgrad-Mdat.Sp2.*sign(Wk2{1});
       else
           GradBlock=BPKF_Extract_Block(WKgrad,[Mdat.Block{:}]);
           %% Initialize with L1 gradients
           GradVal.Wk1=Uncellfun(@(xx)(-Mdat.Sp1.*sign(xx)),Wk1);
           GradVal.Wk2=Uncellfun(@(xx)(-Mdat.Sp2.*sign(xx)),Wk2);
           %% Accumulate through blocks (this can be made more efficient...)
           for iBlock=1:numel(GradBlock)
               pre0=Mdat.PreGroup{iBlock};
               post0=Mdat.PostGroup{iBlock};
               GradVal.Wk1{pre0}=GradVal.Wk1{pre0}+GradBlock{iBlock}*Wk2{post0}';
               GradVal.Wk2{post0}=GradVal.Wk2{post0}+Wk1{pre0}'*GradBlock{iBlock};
           end
       end
    end

%%%%%%%%%%%%%%%%%%%%%%%% Evaluate Regularization in the form:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%    L1coeff*||X-L1mean||_1+(L2coeff/2)*||X-L2mean||^2_2

    L1names=ParamValDep.Reg.L1names;
    if ~isempty(L1names)
    L1coeff=ParamValDep.Reg.L1.coeff;
    L1mean=ParamValDep.Reg.L1.mean;
    for iReg1=1:numel(L1names)
        GradVal.(L1names{iReg1})=GradVal.(L1names{iReg1})-...
            L1coeff{iReg1}.*sign(ParamValDep.ParamStr.(L1names{iReg1})-L1mean{iReg1});
    end
    end

    L2names=ParamValDep.Reg.L2names;
    if ~isempty(L2names)
    L2coeff=ParamValDep.Reg.L2.coeff;
    L2mean=ParamValDep.Reg.L2.mean;
    for iReg2=1:numel(L2names)
        GradVal.(L2names{iReg2})=GradVal.(L2names{iReg2})-...
            L2coeff{iReg2}.*(ParamValDep.ParamStr.(L2names{iReg2})-L2mean{iReg2});
    end
    end
end

if nadamStr.Spec.isClip
    combGrads=0;
    for ii=1:numel(ParamValDep.ParamStr.varOrder)
        combGrads=combGrads+sum(GradVal.(ParamValDep.ParamStr.varOrder{ii}).^2,[1 2]);
    end
    normGrad=sqrt(combGrads);
    %combGrads=cell(1,numel(ParamValDep.ParamStr.varOrder));
%for ii=1:numel(combGrads)
%    combGrads{ii}=GradVal.(ParamValDep.ParamStr.varOrder{ii});
%end
%combGrads=[combGrads{:}];

%    normGrad=sqrt(sum(combGrads.^2,[1 2]));
    if iBatch<=nadamStr.Spec.ClipLength
            nadamStr.normHistory(iBatch)=normGrad;
            clipFactor=1;
    else
%        mean(nadamStr.Spec.clipHistory)
        normThresh=mean(nadamStr.normHistory)*nadamStr.Spec.Clip;
        if normGrad>normThresh
            clipFactor=normThresh/normGrad;
            ParamValDep.censTime(1,iBatch)=normGrad/normThresh;
        else
            clipFactor=1;
        end
        %% Note that history updates use the semi-corrected gradient norm
        if strcmpi(nadamStr.Spec.ClipType(1),'m')
            nadamStr.normHistory=[nadamStr.normHistory(2:end) (clipInnovation*normGrad+(1-clipInnovation)*normThresh)];
        end
    end
else
    clipFactor=1;
end
if ~isempty(ParamValDep.DepStr)
GradVal=BPKF_Dependency_Gradients(ParamValDep,GradVal);
end


mu=nadamStr.Spec.dec1;
v=nadamStr.Spec.dec2;
if nadamStr.Spec.Alg.doRateFunc
    Rate=nadamStr.Spec.Rate(iBatch);
else
Rate=nadamStr.Spec.Rate;
end


Alg=nadamStr.Spec.Alg;
doAR1=Alg.doAR1;
doAR2=Alg.doAR2;


gHat=(1-mu)/(1-mu^iBatch);
mHat=mu/(1-mu^(iBatch+1));
v0=1/(1-v^iBatch);
Rtv0=sqrt(v0);
f1=(Rate*mHat/Rtv0);
f2=(Rate*gHat/Rtv0);

Reg0=nadamStr.Spec.reg/Rtv0;
gradFields=fields(GradVal);
if ~isempty(ParamValDep.DepStr)
depFields=fields(ParamValDep.DepStr.Var);
else
    depFields={};
end
isDep=ismember(gradFields,depFields);


for iVar=1:numel(gradFields)
    pName=gradFields{iVar};
    %isDep=ismember(pName,depFields);
    nTmp=['n',pName];
    mTmp=['m',pName];
    sTmp=[pName,'Sign'];
    
    if isDep(iVar)
    if iscell(GradVal.(pName))
        for iGroup=1:numel(GradVal.(pName))
    gradX=clipFactor*GradVal.(pName){iGroup};
    if doAR1
    nadamStr.(mTmp){iGroup}=mu*nadamStr.(mTmp){iGroup}+(1-mu)*gradX;
    end
    if doAR2
    nadamStr.(nTmp){iGroup}=v*nadamStr.(nTmp){iGroup}+(1-v)*gradX.^2;
    end
    if Alg.doNADAM
        deltaX=((f1*nadamStr.(mTmp){iGroup})+(f2*gradX))./(sqrt(nadamStr.(nTmp){iGroup})+Reg0);
    elseif Alg.doSGD
        deltaX=Rate*gradX;
    elseif Alg.doMom
        deltaX=Rate*nadamStr.(mTmp){iGroup};
    elseif Alg.doRMSprop
        deltaX=Rate.*gradX./sqrt(nadamStr.(nTmp){iGroup}+nadamStr.Spec.reg);
    end
    newVal=ParamValDep.DepStr.Var.(pName){iGroup}+deltaX;
        %% Basis vectors are always positive
        if strcmpi(pName,'Vec')
            ParamValDep.DepStr.Var.(pName){iGroup}=max(0,newVal);%./...
%                max(abs(newVal),[],1);
        else
            tmpSign=ParamValDep.DepStr.Dep.(sTmp){iGroup};
            ParamValDep.DepStr.Var.(pName){iGroup}(:,tmpSign==1)=...
                max(0,newVal(:,tmpSign==1));
            ParamValDep.DepStr.Var.(pName){iGroup}(:,tmpSign==-1)=...
                min(0,newVal(:,tmpSign==-1));
        end
        end
    else
    gradX=clipFactor*GradVal.(pName);
    if doAR1
    nadamStr.(mTmp)=mu*nadamStr.(mTmp)+(1-mu)*gradX;
    end
    if doAR2
    nadamStr.(nTmp)=v*nadamStr.(nTmp)+(1-v)*gradX.^2;
    end
    if Alg.doNADAM
    deltaX=((f1*nadamStr.(mTmp))+(f2*gradX))./(sqrt(nadamStr.(nTmp))+Reg0);
    elseif Alg.doSGD
        deltaX=Rate*gradX;
    elseif Alg.doMom
        deltaX=Rate*nadamStr.(mTmp);
    elseif Alg.doRMSprop
        deltaX=Rate.*gradX./sqrt(nadamStr.(nTmp)+nadamStr.Spec.reg);
    end
    newVal=ParamValDep.DepStr.Var.(pName)+deltaX;

    tmpSign=ParamValDep.DepStr.Dep.(sTmp);
    ParamValDep.DepStr.Var.(pName)(:,tmpSign==1)=...
                max(0,newVal(:,tmpSign==1));
    ParamValDep.DepStr.Var.(pName)(:,tmpSign==-1)=...
                min(0,newVal(:,tmpSign==-1));
    end
    else
        if iscell(GradVal.(pName))
            for iGroup=1:numel(GradVal.(pName))
                    gradX=clipFactor*GradVal.(pName){iGroup};
                    if doAR1
    nadamStr.(mTmp){iGroup}=mu*nadamStr.(mTmp){iGroup}+(1-mu)*gradX;
                    end
                    if doAR2
    nadamStr.(nTmp){iGroup}=v*nadamStr.(nTmp){iGroup}+(1-v)*gradX.^2;
                    end
                    if Alg.doNADAM
                        deltaX=((f1*nadamStr.(mTmp){iGroup})+(f2*gradX))./(sqrt(nadamStr.(nTmp){iGroup})+Reg0);
                    elseif Alg.doSGD
                        deltaX=Rate*gradX;
                    elseif Alg.doMom
                        deltaX=Rate*nadamStr.(mTmp){iGroup};
                    elseif Alg.doRMSprop
                         deltaX=Rate.*gradX./sqrt(nadamStr.(nTmp){iGroup}+nadamStr.Spec.reg);                       
                    end
    ParamValDep.ParamStr.(pName)=ParamValDep.ParamStr.(pName)+deltaX;
            end
        else
    gradX=clipFactor*GradVal.(pName);
            if doAR1
    nadamStr.(mTmp)=mu*nadamStr.(mTmp)+(1-mu)*gradX;
            end
            if doAR2
    nadamStr.(nTmp)=v*nadamStr.(nTmp)+(1-v)*gradX.^2;
            end
            if Alg.doNADAM
                deltaX=((f1*nadamStr.(mTmp))+(f2*gradX))./(sqrt(nadamStr.(nTmp))+Reg0);
            elseif Alg.doSGD
                deltaX=Rate*gradX;
            elseif Alg.doMom
                deltaX=Rate*nadamStr.(mTmp);
            elseif Alg.doRMSprop
                deltaX=Rate.*gradX./sqrt(nadamStr.(nTmp)+nadamStr.Spec.reg);
            end
    ParamValDep.ParamStr.(pName)=ParamValDep.ParamStr.(pName)+deltaX;
        end
    end
end
%if ParamValDep.ParamStr.doReg
%    if isfield(ParamValDep.ParamStr,'Wk1')
%    ParamValDep.ParamStr.W=ParamValDep.ParamStr.W0+ParamValDep.ParamStr.Wk1*ParamValDep.ParamStr.Wk2';
%    end
%end
ParamValDep.ParamStr=BPKF_Apply_Dependency(ParamValDep);
        
        
end