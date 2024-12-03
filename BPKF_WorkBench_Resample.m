function[DIRstr]=BPKF_WorkBench_Resample(DIRstr,Lname,Rname)
saveDir=DIRstr.saveDir;
%% Directory for HCP workbench
WBdir=DIRstr.WorkBench;
%% Original Surface Directory
surfDIR=DIRstr.Surf.Folder;
%% Original Surface
surfDIR_L=DIRstr.Surf.L;
surfDIR_R=DIRstr.Surf.R;
%% Original Sphere
surfDIR_Lsphere=DIRstr.Surf.Lsphere;
surfDIR_Rsphere=DIRstr.Surf.Rsphere;
%% Surface Template Directory
templateDIR=DIRstr.Template.Folder;
%% Surface Template Files
templateDIR_Lsphere=DIRstr.Template.Lsphere;
templateDIR_Rsphere=DIRstr.Template.Rsphere;


BATname='BPKF_Resample_Script.bat';

giiEnd='.surf.gii';


if numel(Lname)<=numel(giiEnd) || ~strcmp(Lname((end-numel(giiEnd)+1):end),giiEnd)
    Lname=strcat(Lname,giiEnd);
end
if numel(Rname)<=numel(giiEnd) || ~strcmp(Rname((end-numel(giiEnd)+1):end),giiEnd)
    Rname=strcat(Rname,giiEnd);
end

sphereLname=['Sphere_',Lname];
sphereRname=['Sphere_',Rname];


surfFiles={surfDIR_L,surfDIR_R,surfDIR_Lsphere,surfDIR_Rsphere};
templateFiles={templateDIR_Lsphere,templateDIR_Rsphere};

for ii=1:numel(surfFiles)
    copyfile(fullfile(surfDIR,surfFiles{ii}),fullfile(saveDir,surfFiles{ii}))
end
for ii=1:numel(templateFiles)
    copyfile(fullfile(templateDIR,templateFiles{ii}),fullfile(saveDir,templateFiles{ii}));
end

fid = fopen(fullfile(saveDir,BATname),'w'); 


CommandCell{1}=['cd ',saveDir];
ReSample_Start=[fullfile(WBdir,'wb_command'),' -surface-resample '];
CommandCell{2}=[ReSample_Start,surfDIR_L,' ',surfDIR_Lsphere,' ',templateDIR_Lsphere, ' BARYCENTRIC ',Lname];
CommandCell{3}=[ReSample_Start,surfDIR_R,' ',surfDIR_Rsphere,' ',templateDIR_Rsphere, ' BARYCENTRIC ',Rname];
CommandCell{4}=[ReSample_Start,surfDIR_Lsphere,' ',surfDIR_Lsphere,' ',templateDIR_Lsphere, ' BARYCENTRIC ',sphereLname];
CommandCell{5}=[ReSample_Start,surfDIR_Rsphere,' ',surfDIR_Rsphere,' ',templateDIR_Rsphere, ' BARYCENTRIC ',sphereRname];


NewNames={Lname,Rname,sphereLname,sphereRname};

fprintf(fid, '%s\n',CommandCell{:}) ;
fclose(fid) ;

for ii=1:numel(NewNames)
fid=fopen(fullfile(saveDir,NewNames{ii}),'w');
fclose(fid);
end

system(fullfile(saveDir,BATname))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Overwrite Surface File Locations to New Values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DIRstr.Surf.Folder=saveDir;
DIRstr.Surf.L=Lname;
DIRstr.Surf.R=Rname;
DIRstr.Surf.Lsphere=sphereLname;
DIRstr.Surf.Rsphere=sphereRname;

end