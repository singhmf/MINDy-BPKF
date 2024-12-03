function[PBTcost,PBTstr]=BPKF_PBT_Evaluate(Out,PBTstr,iBatch)
%% Evaluate cost for population-based training (PBT.cost)
%% PBT.Cost can be a function handle of Out (and iBatch) or a structure specifying weights to be given
%% to kalman predictions ("kal") and free-running predictions ("rec") [scalar or vector of time]

%% REMARK: WILL NEED TO CHANGE THE INDEXING METHOD FOR CV OBJECTIVES

if isa(PBTstr.Cost,"function_handle")
    if nargin(PBTstr.Cost)==1
        PBTcost=PBTstr.Cost(Out);
    elseif nargin(PBTstr.Cost)==2
        PBTcost=PBTstr.Cost(Out,iBatch);
    else
        error('Too many input arguments for PBTstr.Cost')
    end
else
    cTmp=PBTstr.Cost;
    batchInd=unique(max(1,iBatch-(0:(cTmp.Length-1))));
    
    %% batch start/stop used to identify CV indices
        bStart=find(Out.BatchSchedule==min(batchInd));
        bEnd=find(Out.BatchSchedule==iBatch);
        bTmp=Out.BatchSchedule(bStart:bEnd);
        cvMark=abs(bTmp(bTmp<0));

    PBTcost=0;
    costFields={'recE','kalE','recCV','kalCV'};
    IndSet={batchInd,batchInd,cvMark,cvMark};
    for iF=1:numel(costFields)
        cField=costFields{iF};
        if isfield(cTmp,cField)&& ~isempty(cTmp.(cField))
            PBTcost=PBTcost+mean(cTmp.(cField)*mean(Out.(cField)(:,IndSet{iF}),2),'all');
        end
    end
    if isfield(cTmp,'cov')
        PBTcost=PBTcost+cTmp.cov*mean(Out.CovErr(batchInd),'all');
    end
    
end
end
