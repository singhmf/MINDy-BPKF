function[OutMean,OutSTD]=BPKF_Average_Outputs(X,varargin)

tFields={'recCV','kalCV','recE','kalE','recW','recD','recS','recC','recV','covErr','rtQ','Param','StartupTime','runTime'};

if ~iscell(X)
    X=[{X},varargin];
elseif ~isempty(varargin)
    error('Unclear input: either single input of cell or many individual structure inputs');
end
OutMean=X{1};

if nargout==2
    OutSTD=X{1};
end
nX=numel(X);


for iField=1:numel(tFields)
    tmpField=tFields{iField};
    if isfield(OutMean,tmpField)
        cellField=iscell(OutMean.(tmpField));
        for iX=2:nX
            if cellField
                OutMean.(tmpField)=Uncellfun(@plus,OutMean.(tmpField),X{iX}.(tmpField));
            else
                OutMean.(tmpField)=OutMean.(tmpField)+X{iX}.(tmpField);
            end
        end
        if cellField
            OutMean.(tmpField)=Uncellfun(@(xx)(xx/nX),OutMean.(tmpField));
        else
            OutMean.(tmpField)=OutMean.(tmpField)/nX;
        end

        if nargout==2
            if cellField
                OutSTD.(tmpField)=Uncellfun(@(xx,yy)((xx-yy).^2),X{1}.(tmpField),OutMean.(tmpField));
            else
                OutSTD.(tmpField)=(X{1}.(tmpField)-OutMean.(tmpField)).^2;
            end
            for iX=2:nX
                if cellField
                    OutSTD.(tmpField)=Uncellfun(@(xx,yy,zz)(xx+(yy-zz).^2),OutSTD.(tmpField),X{iX}.(tmpField),OutMean.(tmpField));
                else
                    OutSTD.(tmpField)=OutSTD.(tmpField)+(X{iX}.(tmpField)-OutMean.(tmpField)).^2;
                end
            end
            if cellField
                OutSTD.(tmpField)=Uncellfun(@(xx)(sqrt(xx/(nX-1))),OutSTD.(tmpField));
            else
                OutSTD.(tmpField)=sqrt(OutSTD.(tmpField)/(nX-1));
            end
        end


    end
end
end

