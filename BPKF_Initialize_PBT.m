function [PBTstr] = BPKF_Initialize_PBT(PBTstr)
%% Normalize some ease-of-use settings with PBT
disp('Initializing PBT')

VarOrder={'Grad','Param','Kal','Model'};
PBTstr.VarOrder=VarOrder;

makeHoriz=@(xx)(reshape(xx,1,numel(xx)));

if isfield(PBTstr,'lB')&&isstruct(PBTstr.lB)
    lBall=cell(1,numel(VarOrder));
    for ii=1:numel(VarOrder)
        if isfield(PBTstr.lB,VarOrder{ii})
            lBall{ii}=makeHoriz(PBTstr.lB.(VarOrder{ii}));
        end
    end
    PBTstr.lB=[lBall{:}];
end

if isfield(PBTstr,'uB')&&isstruct(PBTstr.uB)
    uBall=cell(1,numel(VarOrder));
    for ii=1:numel(VarOrder)
        if isfield(PBTstr.uB,VarOrder{ii})
            uBall{ii}=makeHoriz(PBTstr.uB.(VarOrder{ii}));
        end
    end
    PBTstr.uB=[uBall{:}];
end
end
