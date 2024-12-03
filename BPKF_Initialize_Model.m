function[ParamValDep,ModelSpec,Hset,Rset,nadamStr]=BPKF_Initialize_Model(ModelSpec,GradSpec,nSet)


%% Optionally shuffle model initialization but use fixed rng seed thereafter
%%      (useful for maintaining a fixed training schedule)
if isfield(ModelSpec,'trainSeed')&&~isempty(ModelSpec.trainSeed)
    isTrainSeed=true;
disp('Shuffling parameter initialization then using fixed seed')
    rng shuffle
else
    isTrainSeed=false;
end


%% 1. Set Defaults
ModelSpec=BPKF_Initialize_Defaults(ModelSpec);

%% 2. Initialization (Error checking) for GradSpec
GradSpec=BPKF_Initialize_GradSpec(GradSpec);

%% 3. Measurement model configuration
[Rset,Hset]=BPKF_Extract_Meas_Model(ModelSpec,nSet);

%% 4. Parameter Initialization
[ParamStr,nadamStr]=BPKF_Initialize_Param(ModelSpec);
if ParamStr.isFixed.Qchol
ParamStr.cholQmat=ParamStr.cholQorig;%chol(ParamStr.Q,'lower');
else
    ParamStr.cholQmat=diag(ParamStr.Qchol);
end

[DepStr,ParamStr,nadamStr]=BPKF_Initialize_Dependency(ModelSpec,ParamStr,nadamStr);
if ~isempty(DepStr)
ParamStr=BPKF_Apply_Dependency(ParamStr,DepStr);
end
nadamStr.Spec=GradSpec;

if isfield(ModelSpec,'Start')
    startVars=fields(ModelSpec.Start);
    disp(['Using pre-specified initialization for: ',[startVars{:}]])
    for ii=1:numel(startVars)
        ParamStr.(startVars{ii})=ModelSpec.Start.(startVars{ii});
    end
end
if GradSpec.isClip
    nadamStr.normHistory=zeros(1,GradSpec.ClipLength);
end

[ParamStr,Reg,nadamStr]=BPKF_Initialize_Reg(ModelSpec,ParamStr,nadamStr);
ParamValDep.Reg=Reg;
ParamValDep.ParamStr=ParamStr;ParamValDep.DepStr=DepStr;
%% Return to fixed training schedule (rng seed) thereafter
if isTrainSeed
    rng(ModelSpec.trainSeed);
end
end