function[GradSpec]=BPKF_Initialize_GradSpec(GradSpec)
%% Initialization and error-checking for the gradient algorithm specification (GradSpec)

%% Gradient clipping threshold is relative the average gradient norm over [k] initial conditions
GradSpec.isClip=isfield(GradSpec,'Clip')&& isnumeric(GradSpec.Clip);
if GradSpec.isClip
    if ~isfield(GradSpec,'ClipLength')
    error('To use gradient clipping, the duration to use for establishing the clip threshold must be set')
    end
    if ~isfield(GradSpec,'GradStart')
    GradSpec.GradStart=0;
    end
    if ~isfield(GradSpec,'ClipType')||~ismember(upper(GradSpec.ClipType(1)),'MS')
        disp('No ClipType specified (either start or moving). Assuming clip based upon start (initial norms)')
        GradSpec.ClipType='Start';
    end
    GradSpec.clipHistory=zeros(1,GradSpec.ClipLength);
end


if isfield(GradSpec,'Type')
    if isfield(GradSpec,'type')
        error('Use either "type" or "Type" for GradSpec not both')
    else
        disp('Replacing "Type" with "type" in GradSpec')
        GradSpec.type=GradSpec.Type;
        GradSpec=rmfield(GradSpec,'Type');
    end
end

gradTypeNames={'NADAM','ADAM','RMSprop','SGD','Mom','NAG','AdaMax'};
for ii=1:numel(gradTypeNames)
    GradSpec.Alg.(strcat('do',gradTypeNames{ii}))=false;
end


if ~isfield(GradSpec,'type')||isempty(GradSpec.type)
    disp('Assuming NADAM updates')
    GradSpec.Alg.doNADAM=true;
    GradSpec.Alg.doAR1=true;
    GradSpec.Alg.doAR2=true;
elseif strcmpi(GradSpec.type,'nadam')
    GradSpec.Alg.doNADAM=true;
    GradSpec.Alg.doAR1=true;
    GradSpec.Alg.doAR2=true;
elseif strcmpi(GradSpec.type,'adam')
    GradSpec.Alg.doADAM=true;
    GradSpec.Alg.doAR1=true;
    GradSpec.Alg.doAR2=true;
elseif strcmpi(GradSpec.type,'sgd')
    GradSpec.Alg.doSGD=true;
    GradSpec.Alg.doAR1=false;
    GradSpec.Alg.doAR2=false;
elseif strcmpi(GradSpec.type,'RMSprop')
    GradSpec.Alg.doRMSprop=true;
    GradSpec.Alg.doAR1=false;
    GradSpec.Alg.doAR2=true;
elseif strcmpi(GradSpec.type(1:3),'mom')
    GradSpec.Alg.doMom=true;
    GradSpec.Alg.doAR1=true;
    GradSpec.Alg.doAR2=false;
elseif strcmpi(GradSpec.type,'nag')
    GradSpec.Alg.doNAG=true;
    GradSpec.Alg.doAR1=false;
    GradSpec.Alg.doAR2=false;
elseif strcmpi(GradSpec.type,'AdaMax')
    GradSpec.Alg.doAdaMax=true;
    GradSpec.Alg.doAR1=true;
    GradSpec.Alg.doAR2=false;
else
    disp('Supported Gradient Algorithms:')
    disp(gradTypeNames)
    error('Unrecognized gradient algorithm')
end
if isa(GradSpec.Rate,'function_handle')
    GradSpec.Alg.doRateFunc=true;
else
    GradSpec.Alg.doRateFunc=false;
end
end


