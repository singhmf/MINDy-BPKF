function[PBTstrCell]=BPKF_PBT_Initialize_Group(tmpFolder,nWorker,PBTsetup)

groupInd=randi(999999,1);
mkdir(fullfile(tmpFolder),strcat('Worker_Group_',num2str(groupInd)));

PBTstrCell=repmat({PBTsetup},1,nWorker);
for ii=1:nWorker
    PBTstrCell.Worker=ii;
    PBTstrCell.Group=groupInd;
    PBTstrCell.isLeader=ii==1;
end
end

