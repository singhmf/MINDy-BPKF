function[Dat]=BPKF_Neutral2cell(Dat,varargin)
%% If Dat is a cell array, take the varargin indices of it
if iscell(Dat)
    Dat=Dat{varargin{:}};
end
end