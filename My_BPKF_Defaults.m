




tmptmpDir=strsplit(cd,'\');tmptmpDir=tmptmpDir{3};

DIRstr.WorkBench=fullfile('C:\Users\',tmptmpDir,'\Desktop\HCP\workbench\bin_windows64');
DIRstr.FieldTrip=fullfile('C:\Users\',tmptmpDir,'\Documents\fieldtrip-master\fieldtrip-master\');
DIRstr.MRI=             fullfile('C:\Users\',tmptmpDir,'\Desktop\Stim Data\9988\freesurfer\sub-9988\mri\');
DIRstr.saveDir=         fullfile('C:\Users\',tmptmpDir,'\Desktop\Stim Data\Tmp_Test_Resample_Script\');
DIRstr.Surf.Folder=     fullfile('C:\Users\',tmptmpDir,'\Desktop\Stim Data\Resampling\');
DIRstr.Template.Folder=     fullfile('C:\Users\',tmptmpDir,'\Desktop\Stim Data\Resampling\');



DIRstr.Template.Lsphere=    'fsaverage5_lh.sphere.surf.gii';
DIRstr.Template.Rsphere=    'fsaverage5_rh.sphere.surf.gii';
DIRstr.Surf.L=          'sub-9988-L_midthickness.surf.gii';
DIRstr.Surf.R=          'sub-9988-R_midthickness.surf.gii';
DIRstr.Surf.Lsphere=    '9988_lh.sphere.surf.gii';
DIRstr.Surf.Rsphere=    '9988_rh.sphere.surf.gii';

Fid_Nz=[-5.5 105.5 -2.5];
Fid_LPA=[-74.5 -4.5 -18.5];
Fid_RPA=[78.5 11.5 -20.5];

clear tmptmpDir