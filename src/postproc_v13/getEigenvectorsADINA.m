function [eVec] = getEigenvectorsADINA(porFN,nNode,startMode,endMode,param)

% Read eigenvectors from adina porthole file
% (endMode-startMode+1) eigenvectors will be obtained
% Input:
%        porFN     : the porthole file name
%        nNode     : the number of nodes
%        param     : struct containing parameters 
%                    such as startMode, endMode, ADINA path, ...
% Output:
%        eVec      : nNode*3(x,y,z)*nMode

% remove .por extension
bodyFN = porFN(1:length(porFN)-4);

% write .plo file
ploFN = strcat(bodyFN,'_evec.plo');
listFN = strcat(bodyFN,'_evec.list');

fid = fopen(ploFN,'w');
fprintf(fid,'%s\n','LOADPORTHOLE OPERATIO=CREATE FILE=,');
fprintf(fid,'''');
fprintf(fid,'%s',porFN);
fprintf(fid,''',\n');
fprintf(fid,'%s\n','     TAPERECO=0 DUMPFORM=NO PRESCAN=NO RANGE=ALL TIMESTAR=0.00000000000000,');
fprintf(fid,'%s\n','     TIMEEND=0.00000000000000 STEPSTAR=0 STEPEND=0 STEPINCR=1,');
fprintf(fid,'%s\n','     ZOOM-MOD=0 INITIAL-=AUTOMATIC CPSTART=1 CPEND=0');
fprintf(fid,'*\n');
fprintf(fid,'%s','RESPRANGE MODE-SHAPE NAME=DEFAULT MODESTAR=');
fprintf(fid,'%s',num2str(startMode));
fprintf(fid,'%s',' MODEEND=');
fprintf(fid,'%s',num2str(endMode));
fprintf(fid,',\n');
fprintf(fid,'%s\n','     REFTIME=LATEST');
fprintf(fid,'*\n');
fprintf(fid,'%s\n','FILELIST OPTION=FILE FILE=,');
fprintf(fid,'''');
fprintf(fid,'%s',listFN);
fprintf(fid,'''\n');
fprintf(fid,'*\n');
fprintf(fid,'%s\n','ZONELIST ZONENAME=WHOLE_MODEL RESULTGR=DEFAULT SMOOTHIN=DEFAULT,');
fprintf(fid,'%s\n','     RESULTCO=DEFAULT RESPOPTI=RESPRANGE RESPONSE=DEFAULT,');
fprintf(fid,'%s\n','     RESPRANG=DEFAULT VARIABLE=X-EIGENVECTOR Y-EIGENVECTOR,');
fprintf(fid,'%s\n','     Z-EIGENVECTOR,');
fprintf(fid,'%s\n','     X-EIGENVECTOR_ROTATION Y-EIGENVECTOR_ROTATION Z-EIGENVECTOR_ROTATION');
fprintf(fid,'%s\n','****');
fprintf(fid,'%s\n','*** EXIT SAVE=UNKNOWN PROMPT=NO IMMEDIAT=NO');
fclose(fid);

% run AUI to get eigenvectors
runAUI = sprintf('%s %s %s',param.auiEXE,param.auiOPTION,ploFN);
system(runAUI);

% read eigenvectors
fid = fopen(listFN,'r');
eVec = zeros(nNode,6,endMode-startMode+1);
for i=1:endMode-startMode+1
    while 1
        ss = fgetl(fid);
        if strfind(ss,'Mode_number')
            break;
        end
    end
    ss = fgetl(fid);
    for j=1:nNode
        ss = fgetl(fid);
        raw = textscan(ss,'%*s%*d%f%f%f%f%f%f');
        eVec(j,1,i) = raw{1};
        eVec(j,2,i) = raw{2};
        eVec(j,3,i) = raw{3};
        eVec(j,4,i) = raw{4};
        eVec(j,5,i) = raw{5};
        eVec(j,6,i) = raw{6};
    end
end
fclose(fid);

