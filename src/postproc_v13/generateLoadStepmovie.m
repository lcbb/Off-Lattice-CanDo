function [] = generateLoadStepmovie(tarDIR,bodyFN,param,matPATH,outPATH)

movieDIR = 'movieLoadstep';
load(matPATH,'FEModel');

node = FEModel.node(:,1:3);
if(isfield(FEModel,'rigidLinks'))
    elem = [ ...
            FEModel.beamsNormal; ...
            FEModel.beamsNick; ...
            FEModel.rigidLinks; ...
           ];
else
    elem = [ ... 
            FEModel.beamsNormal; ...
            FEModel.beamsNick; ...
           ];
end
elemBeam = [ ... 
            FEModel.beamsNormal; ...
            FEModel.beamsNick; ...
           ];   
elemRender = [ ...
               FEModel.beamsNormal; ...
               FEModel.beamsNick; ...
%                FEModel.shrinkBeamsNormal; ...
%                FEModel.shrinkBeamsNick; ...
             ];
if(isfield(FEModel,'rigidLinks'))
    elemRenderTruss = [ ...
                       FEModel.rigidLinks; ...
                       FEModel.entSprings; ...
                      ];
else
    elemRenderTruss = [FEModel.entSprings];
end
nNode = size(node,1);
nElem = size(elem,1);
nElemBeam = size(elemBeam,1);
nElemRender = size(elemRender,1);
nElemRenderTruss = size(elemRenderTruss,1);
nAlignHJE = size(FEModel.alignHJE,1);
nAlignDSDNA = size(FEModel.alignDSDNA,1);
nTruss = size(FEModel.entSprings,1);

matFN = fullfile(tarDIR,strcat(bodyFN,'_LoadStepDisplVel.mat'));
% if exist(matFN,'file')
%     load(matFN);
% else
    % read displacements from ADINA .por file
    fid = fopen(fullfile(tarDIR,strcat(bodyFN,'_getLoadStepDispl.plo')),'w');
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
    fprintf(fid,'%s',fullfile(tarDIR,strcat(bodyFN,'_loadstep_displVel.txt')));
    fprintf(fid,'''\n');
    fprintf(fid,'RESPRANGE LOAD-STEP NAME=DEFAULT TSTART=EARLIEST TEND=LATEST,\n');
    fprintf(fid,'     INCREMEN=AVAILABLE INTERPOL=NO NSKIP=0\n');
    fprintf(fid,'ZONELIST ZONENAME=WHOLE_MODEL RESULTGR=DEFAULT SMOOTHIN=DEFAULT,\n');
    fprintf(fid,' RESULTCO=DEFAULT RESPOPTI=RESPRANGE RESPONSE=DEFAULT,\n');
    fprintf(fid,' RESPRANG=DEFAULT VARIABLE=X-DISPLACEMENT Y-DISPLACEMENT,\n');
    fprintf(fid,' Z-DISPLACEMENT X-VELOCITY Y-VELOCITY Z-VELOCITY\n');
    fprintf(fid,'EXIT SAVE=NO IMMEDIATE=YES\n');
    fclose(fid);
    runAUI = sprintf('%s %s %s',param.auiEXE,param.auiOPTION,fullfile(tarDIR,strcat(bodyFN,'_getLoadStepDispl.plo')));
    system(runAUI);

    % read displacements (global Coordinates)
    fid = fopen(fullfile(tarDIR,strcat(bodyFN,'_loadstep_displVel.txt')),'r');
    displ_vel = zeros(nNode,6,1);
    iTime = 0;
    while 1
        ss = fgetl(fid);
        if ~ischar(ss)
            break;
        else
            if strfind(ss,'Time ')
                iTime = iTime + 1;
                [solutionTime(iTime,1)]=strread(ss,'%*s%f');
                ss = fgetl(fid); % skip blank line
                for kkk=1:nNode
                    ss = fgetl(fid);
                    [u1,u2,u3,u4,u5,u6] = strread(ss,'%*s%*s%s%s%s%s%s%s');
                    if isempty(strfind(char(u1),'E'))
                        displ_vel(kkk,1,iTime) = 0.0;
                    else
                        displ_vel(kkk,1,iTime) = str2double(char(u1));
                    end
                    if isempty(strfind(char(u2),'E'))
                        displ_vel(kkk,2,iTime) = 0.0;
                    else
                        displ_vel(kkk,2,iTime) = str2double(char(u2));
                    end
                    if isempty(strfind(char(u3),'E'))
                        displ_vel(kkk,3,iTime) = 0.0;
                    else
                        displ_vel(kkk,3,iTime) = str2double(char(u3));
                    end
                    if isempty(strfind(char(u4),'E'))
                        displ_vel(kkk,4,iTime) = 0.0;
                    else
                        displ_vel(kkk,4,iTime) = str2double(char(u4));
                    end
                    if isempty(strfind(char(u5),'E'))
                        displ_vel(kkk,5,iTime) = 0.0;
                    else
                        displ_vel(kkk,5,iTime) = str2double(char(u5));
                    end
                    if isempty(strfind(char(u6),'E'))
                        displ_vel(kkk,6,iTime) = 0.0;
                    else
                        displ_vel(kkk,6,iTime) = str2double(char(u6));
                    end
                end
            end
        end
    end
    fclose(fid);
    save(matFN,'displ_vel','solutionTime');
% end

% Read displacements
displ_vel = displ_vel(:,:,1:param.skiploadstep:end);
nSnapshot = size(displ_vel,3);
for i=1:param.framerate
    displ_vel(:,:,nSnapshot+i) = displ_vel(:,:,nSnapshot);
end
nSnapshot = size(displ_vel,3);

% Project onto the principal axes of the converged solution
dNode = node + displ_vel(:,1:3,end);
[pAxes,pNode] = pca(dNode);
pNodeSnap = zeros(nNode,3,nSnapshot);
pNode = node*pAxes;
for i=1:nSnapshot
    pNodeSnap(:,:,i) = pNode + displ_vel(:,1:3,i)*pAxes;
end
if dot(cross(pAxes(:,1),pAxes(:,2)),pAxes(:,3))<0
    pNode(:,3) = -pNode(:,3);
    pNodeSnap(:,3,:) = -pNodeSnap(:,3,:);
end

% Read strain energy
matFN = fullfile(tarDIR,strcat(bodyFN,'_LoadStepStrE.mat'));
nelem = size(elemBeam,1);
%solutionTime = solutionTime(1:param.skiploadstep:end);
[strE,trussSS,trussEE,misalignHJE,misalignDSDNA,reaction,solutionTime_outfile] = ...
    procAdinaOut(outPATH,nelem,nTruss,nAlignHJE,nAlignDSDNA,param.timeStep);
assert(max(abs((solutionTime(2:end)-solutionTime_outfile)./solutionTime_outfile)) < 1e-4);
strE = strE / param.KbT;
strE = cat(2,zeros(nelem,1),strE);
trussStrE = getTrussStrE(matPATH,param,trussSS,trussEE);
trussStrE = trussStrE / param.KbT;
trussStrE = cat(2,zeros(nTruss,1),trussStrE);
save(matFN,'strE','trussStrE','misalignHJE','misalignDSDNA','reaction','solutionTime_outfile');
strE = cat(2,strE,repmat(strE(:,end),1,nSnapshot-size(strE,2)));

% Read triads
triad = cell(nNode,nSnapshot);
triad(:,1) = getTriad(matPATH, fullfile(tarDIR,strcat(bodyFN,'.por')), []);
triad(:,2:numel(solutionTime)) = getTriad(matPATH, fullfile(tarDIR,strcat(bodyFN,'.por')), solutionTime(2:end));
triad(:,numel(solutionTime)+1:end) = repmat(triad(:,numel(solutionTime)), 1, nSnapshot-numel(solutionTime));
for loop = 1:nSnapshot
    for i = 1:nNode
        triad{i,loop} = pAxes' * triad{i,loop};
    end
end
if(dot(cross(pAxes(:,1),pAxes(:,2)),pAxes(:,3))<0)
    for loop = 1:nSnapshot
        for i = 1:nNode
            triad{i,loop}(3,:) = -triad{i,loop}(3,:);
        end
    end
end

% Data array ---> BILD & PNG
solutionShape2BildPng(matPATH,pNodeSnap,'_loadstep',[],strE,param.strE_range,triad,'jet',param,movieDIR);

% PNG ---> AVI
[~,prefixFN] = fileparts(matPATH);
prefixFN = prefixFN(1:end-8);   % remove the suffix '_FEModel'
assert(~isempty(prefixFN));
png2avi(fullfile(tarDIR,movieDIR,strcat(prefixFN,'_loadstep')), size(pNodeSnap,3), param.framerate, param);

end