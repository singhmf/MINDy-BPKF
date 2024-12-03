function[PBTstr]=BPKF_PBT_Write(PBTcost,ParamValDep,nadamStr,ModelSpec,KalSpec,PBTstr,iBatch)

    workerInd=PBTstr.Worker;
    TmpDir=PBTstr.Dir;
    iRep=ceil(iBatch/PBTstr.SetLength);
    valEval=PBTcost;
    if ~exist(TmpDir,'dir')
        %% Pause so not all workers try to write at once
        pause(randi(1000)/5000)
        if ~exist(TmpDir,'dir')
            mkdir(TmpDir);
        end
    end


    if isfield(PBTstr,'Group')&&~isempty(PBTstr.Group)
    fName=fullfile(TmpDir,num2str(PBTstr.Group),strcat(num2str(iRep),'_',num2str(valEval,8),'_',num2str(workerInd),'.mat'));
    else
    fName=fullfile(TmpDir,strcat(num2str(iRep),'_',num2str(valEval,8),'_',num2str(workerInd),'.mat'));
    end
    PBTstr0=PBTstr;
   %% Grad, Param, Kal, Model
   HyperVec=cell(1,4);
   HyperInds=cell(1,4);
   CellInds=cell(1,4);
   
   tmp_notUsed=@(xx,yy)(~isfield(xx,yy)||isempty(xx.(yy)));
   

   
   if ~tmp_notUsed(PBTstr.Var,'Grad')
   [HyperVec{1},HyperInds{1},CellInds{1}]=BPKF_PBT_Vectorize(nadamStr,PBTstr.Var.Grad);
   end
   if ~tmp_notUsed(PBTstr.Var,'Param')
   [HyperVec{2},HyperInds{2},CellInds{2}]=BPKF_PBT_Vectorize(ParamValDep,PBTstr.Var.Param);
   end
   if ~tmp_notUsed(PBTstr.Var,'Kal')
   [HyperVec{3},HyperInds{3},CellInds{3}]=BPKF_PBT_Vectorize(KalSpec,PBTstr.Var.Kal);
   end
   if ~tmp_notUsed(PBTstr.Var,'Model')
   [HyperVec{4},HyperInds{4},CellInds{4}]=BPKF_PBT_Vectorize(ModelSpec,PBTstr.Var.Model);
   end

   HyperVec=[HyperVec{:}];


   %% First rep: cost is abs error; other reps: delta error
   if iRep==1 %#ok<IFBDUP> 
       PB2err=PBTcost;
   else
       PB2err=PBTcost;%-PBTstr.Hist.postCost(iRep-1);%-PBTcost;
   end
   
    if isfield(PBTstr,'pendingFile')&&~isempty(PBTstr.pendingFile)
        delete(PBTstr.pendingFile)
        PBTstr.pendingFile=[];
    end    
    save(fName,'ParamValDep','PBTstr0','nadamStr','ModelSpec','HyperVec','HyperInds','CellInds','PB2err');
    PBTstr.HyperInds=HyperInds;
    PBTstr.CellInds=CellInds;

end
