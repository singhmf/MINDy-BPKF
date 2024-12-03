function[ModelSpec]=BPKF_Initialize_Defaults(ModelSpec)


PairSets={  {'ErrMat','L2'},...
            {'RandScale',10^-4},...
            {'minS',.25},...
            {'maxS',2},...
            {'minD',.4},...
            {'maxD',.8}};

for ii=1:numel(PairSets)
    if (~isfield(ModelSpec,PairSets{ii}{1}))||isempty(ModelSpec.(PairSets{ii}{1}))
        ModelSpec.(PairSets{ii}{1})=PairSets{ii}{2};
        disp(['Setting ',PairSets{ii}{1},' equal to ',num2str(PairSets{ii}{2})])
    end
end
end