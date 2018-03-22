% This function is modified from
% Ln 571-610
%
% Two new files will be created in the working directory:
%     - *_getDispl.plo
%     - *_displ.txt
%
function displ = por2displ(porPATH,matPATH,param)

% Read in the finite element model
load(matPATH,'FEModel');
node = FEModel.node;
clear('FEModel');

% Get the prefix of the filenames
[tarDIR,bodyFN] = fileparts(porPATH);

% Generate and run the ADINA *.plo file
% Result: a *.txt file is created
fid = fopen(fullfile(tarDIR,strcat(bodyFN,'_getDispl.plo')),'w');
fprintf(fid,'PLSYSTEM SYSTEM=OPENGL\n');
fprintf(fid,'LOADPORTHOLE OPERATIO=CREATE FILE=,\n');
fprintf(fid,'''');
fprintf(fid,'%s',fullfile(tarDIR,strcat(bodyFN,'.por')));
fprintf(fid,''',\n');
fprintf(fid,'%s\n','     TAPERECO=0 DUMPFORM=NO PRESCAN=NO RANGE=ALL TIMESTAR=0.00000000000000,');
fprintf(fid,'%s\n','     TIMEEND=0.00000000000000 STEPSTAR=0 STEPEND=0 STEPINCR=1,');
fprintf(fid,'%s\n','     ZOOM-MOD=0 INITIAL-=AUTOMATIC CPSTART=1 CPEND=0');
fprintf(fid,'%s\n','FILELIST OPTION=FILE FILE=,');
fprintf(fid,'''');
fprintf(fid,'%s',fullfile(tarDIR,strcat(bodyFN,'_displ.txt')));
fprintf(fid,'''\n');
fprintf(fid,'ZONELIST ZONENAME=WHOLE_MODEL RESULTGR=DEFAULT SMOOTHIN=DEFAULT,\n');
fprintf(fid,' RESULTCO=DEFAULT RESPOPTI=RESPONSE RESPONSE=DEFAULT,\n');
fprintf(fid,' RESPRANG=DEFAULT VARIABLE=X-DISPLACEMENT Y-DISPLACEMENT,\n');
fprintf(fid,' Z-DISPLACEMENT X-ROTATION Y-ROTATION Z-ROTATION\n');
fprintf(fid,'EXIT SAVE=NO IMMEDIATE=YES\n');
fclose(fid);
runAUI = sprintf('%s %s %s',param.auiEXE,param.auiOPTION,fullfile(tarDIR,strcat(bodyFN,'_getDispl.plo')));
system(runAUI);

% Read displacements (in global coordinates) from the newly created *.txt
% file
fid = fopen(fullfile(tarDIR,strcat(bodyFN,'_displ.txt')),'r');
displ = zeros(size(node,1),6);
while 1
    ss = fgetl(fid);
    if strfind(ss,'Time ')
        break;
    end
end
ss = fgetl(fid);
for kkk=1:size(node,1)
    ss = fgetl(fid);
    [displ(kkk,1),displ(kkk,2),displ(kkk,3),displ(kkk,4),displ(kkk,5),displ(kkk,6)]=strread(ss,'%*s%*d%f%f%f%f%f%f');
end
fclose(fid);

end