function[X]=BPKF_Convert_Cell_Size(X,nRep,varName)
%% Reformat variable into nRep repeated cells (if not already)
if iscell(X)
    if numel(X)==1
        X=repmat(X,1,nRep);
    elseif numel(X)~=nRep
        error(['Field ',varName,' has ',num2str(numel(X)),' elements instead of ',num2str(nRep)]);
    end
else
    X=repmat({X},1,nRep);
end
