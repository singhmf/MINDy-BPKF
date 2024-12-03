function[Rset,Hset]=BPKF_Extract_Meas_Model(ModelSpec,nSet)
%% Reformat and retrieve parameters
Rset=ModelSpec.R;
Hset=ModelSpec.H;
if ~iscell(Hset)
    Hset=repmat({Hset},1,nSet);
elseif numel(Hset)==1
    Hset=repmat(Hset,1,nSet);
end
if ~iscell(Rset)
    Rset=repmat({Rset},1,nSet);
elseif numel(Rset)==1
    Rset=repmat(Rset,1,nSet);
end


if isfield(ModelSpec,'Hkron')
Hset=Uncellfun(@(xx)(kron(ModelSpec.Hkron,xx)),Hset);
end

end