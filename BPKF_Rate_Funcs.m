function[Out]=BPKF_Rate_Funcs(LearnType,maxBatch,baseRate,varargin)


revVar=1-varargin{1};
if strcmpi(LearnType(1:3),'Lin')
    disp('Interpreting hyperparameter as minimum proportion')
    Out=@(xx)(baseRate*(revVar*(1-xx./maxBatch)+varargin{1}));
elseif strcmpi(LearnType(1:3),'cos')
    disp('Interpreting hyperparameter as minimum proportion')
    Out=@(xx)(baseRate*(varargin{1}+(revVar/2)*(1+cos(xx*pi/maxBatch))));
elseif strcmpi(LearnType(1:3),'exp')
    disp('Interpreting hyperparameter as exp decay rate')
    Out=@(xx)(baseRate*(exp(-varargin{1}*xx/maxBatch)));
elseif strcmpi(LearnType(1:3),'cyc')
    disp('Interpreting hyperparameters as minimum prop., period')
    Out=@(xx)(baseRate*(varargin{1}+revVar*(1-2*abs(round(xx./varargin{2})-xx./varargin{2}))));
elseif strcmpi(LearnType(1:3),'rex')
    disp('Interpreting hyperparameter as minimum prop.')
    Out=@(xx)(baseRate*(varargin{1}+revVar*2*(1-xx./maxBatch)./(2-xx./maxBatch)));
else
    error('Unrecognized rate function specified')
end
end