function[measComb,guessComb,ModelComb,ParStr]=BPKF_Setup_CV_Inputs(measTrain,guessTrain,ModelSpecTrain,measCV,guessCV,ModelSpecCV,varargin)
%% Combines training and CV data and measurement models (all else is assumed the same) 
%% Remark: Only gets the H. R subfields from ModelSpec CV

if ~isempty(varargin)
    ParStr=varargin{1};
end

measComb=[tmpVec(measTrain),tmpVec(measCV)];
guessComb=[tmpVec(guessTrain),tmpVec(guessCV)];

ModelComb=ModelSpecTrain;

Rtrain=tmpVec(ModelSpecTrain.R);
Rcv=tmpVec(ModelSpecCV.R);

Htrain=tmpVec(ModelSpecTrain.H);
Hcv=tmpVec(ModelSpecCV.H);

if numel(Rtrain)~=numel(conv_cell(measTrain))
Rtrain=repmat(Rtrain,1,numel(measTrain));
end
if numel(Rcv)~=numel(conv_cell(measCV))
Rtrain=repmat(Rcv,1,numel(measCV));
end
if numel(Htrain)~=numel(conv_cell(measTrain))
Htrain=repmat(Htrain,1,numel(measTrain));
end
if numel(Hcv)~=numel(conv_cell(measCV))
Hcv=repmat(Hcv,1,numel(measCV));
end

ModelComb.R=[Rtrain,Rcv];
ModelComb.H=[Htrain,Hcv];
ParStr.CV.mark=numel(conv_cell(measTrain))+(1:numel(conv_cell(measCV)));

end

    function[tmpOut]=tmpVec(tmpIn)
        if ~iscell(tmpIn)
tmpOut={tmpIn};
        else
tmpOut=reshape(tmpIn,1,numel(tmpIn));
        end
    end

    function[tmpIn2]=conv_cell(tmpIn2)
    if ~iscell(tmpIn2)
        tmpIn2={tmpIn2};
    end
    end