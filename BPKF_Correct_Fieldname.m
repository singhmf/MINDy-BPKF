function[oStr]=BPKF_Correct_Fieldname(oStr,oField,varargin)
%% This function searches/coorects for errors of capitalization in specifying field oField in structure oStr
%% Also applies to mis-spellings in varargin including capitalization errors


foundMatch=0;

tmpFields=fields(oStr);
if ~iscell(tmpFields)
    tmpFields={tmpFields};
end
matchField=and(strcmpi(tmpFields,oField),~strcmp(tmpFields,oField));
if sum(matchField)>1
    error(['Overlapping name choices (+/- capitalization) for field ',oField])
elseif sum(matchField)==1
    disp(['Converting fieldname',tmpFields{matchField},' to ',oField]);
    oStr.(oField)=oStr.(tmpFields{matchField});
    oStr=rmfield(oStr,tmpFields{matchField});
    foundMatch=foundMatch+1;
end

tmpFields=fields(oStr);
if ~iscell(tmpFields)
    tmpFields={tmpFields};
end

if ~isempty(varargin)
    for ii=1:numel(varargin)
        matchField=and(strcmpi(tmpFields,varargin{ii}),~strcmp(tmpFields,oField));
        foundMatch=sum(matchField)+foundMatch;
            if sum(matchField)>=1 && foundMatch>1
                error(['Should not be supplying fieldname ',tmpFields{find(matchField,1)},'; instead use: ',oField]);
            else
                tmpName=tmpFields{matchField};
            oStr.(oField)=oStr.(tmpName);
            disp(['Converting fieldname ',tmpName,' to ',oField])
            foundMatch=foundMatch+1;
            oStr=rmfield(oStr,tmpName);
            end
    end
end
end