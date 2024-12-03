function[Dat,leadfield,elec]=BPKF_Mean_Reference_Data(Dat,leadfield,elec,refName)
%% Takes data referenced to channel refName and returns data referenced to mean electric potential
%% Also reduces/reorders LF based upon channels in elec
%% Leadfield MUST contain the contribution of reference channel.
%% Do this BEFORE computing R--so that the channel ordering of R will match new reference scheme
%%      (and include the reference channel)--will have to later mean-recenter R (this should also be implicit in dim. red. of L)


reOrder_elec=zeros(numel(elec.label)+1,1);

for ii=1:numel(elec.label)
    reOrder_elec(ii)=find(strcmpi(elec.label{ii},leadfield.label));
end





reOrder_elec(end)=find(strcmpi(refName,leadfield.label));

nonEmpty=~cellfun(@isempty,leadfield.leadfield);

for ii=find(nonEmpty)
    tmp=leadfield.leadfield{ii}(reOrder_elec,:);
    leadfield.leadfield{ii}=tmp-mean(tmp,1);
end

if ~iscell(Dat)
    Dat=[Dat;zeros(1,size(Dat,2))];
    Dat=Dat-mean(Dat,1);
else
for ii=1:numel(Dat)
    Dat{ii}=[Dat{ii};zeros(1,size(Dat{ii},2))];
    Dat{ii}=Dat{ii}-mean(Dat{ii},1);
end
end
elec.label{end+1}=refName;
end







