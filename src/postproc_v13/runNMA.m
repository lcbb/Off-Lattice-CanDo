% Step 3. Get the properties for the nodes and elements
% Step 3.1 Get the RMSF value, eigenvalue (amplitude), and eigenvector 
%          (shape) for each node
%
function [eVal,eVec,nodalRMSF] = runNMA(porPATH,matPATH,param)

% Get the prefix of the filenames
[tarDIR,bodyFN] = fileparts(porPATH);
initFNbodyNLSA = bodyFN;
if(numel(bodyFN)>5 && strcmp(bodyFN(end-4:end),'_NLSA'))
    bodyFN = bodyFN(1:end-5);
end
initFNbodyNMA = strcat(bodyFN,'_NMA');

% Load topology of the FE model
load(matPATH,'FEModel');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform NMA
% Modified from runCanDo.m, Ln 114-174
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create *.in
generateAdinaIn4NMA(tarDIR, initFNbodyNLSA, initFNbodyNMA, param, FEModel);

% *.in ---> *.dat
adinaIN = fullfile(tarDIR, strcat(initFNbodyNMA,'.in'));
runAUI = sprintf('%s %s %s',param.auiEXE,param.auiOPTION,adinaIN);
system(runAUI);

% *.dat ---> *.out & *.por
adinaDAT = fullfile(tarDIR, strcat(initFNbodyNMA,'.dat'));
runADINA = sprintf('%s %s %s',param.adinaEXE,param.adinaOPTION,adinaDAT);
system(runADINA);

% Check convergence
adinaMSG = fullfile(tarDIR, strcat(initFNbodyNMA,'.msg'));
fid = fopen(adinaMSG);
ss = fgetl(fid);
while(ischar(ss))
    if(strfind(ss,'ENDCODE=1'))
        error('CanDo analysis was partly unsuccessful. The predicted folded shape may be unstable.');
    end
    ss = fgetl(fid);
end
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get eigenvalues and eigenvectors
% Modified from renderNMABildFormat.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Maximum number of normal modes
nNode = size(FEModel.node,1);
startMode = param.startMode;
endMode = min(param.endMode,6*nNode);

% Get eigenvalues
adinaOUT = fullfile(tarDIR,strcat(initFNbodyNMA,'.out'));
eVal = getEigenvaluesADINA(adinaOUT,startMode,endMode);

% Get eigenvectors
adinaPOR = fullfile(tarDIR,strcat(initFNbodyNMA,'.por'));
eVec = getEigenvectorsADINA(adinaPOR,nNode,startMode,endMode,param);

% calculate nodal RMSF
MSF = zeros(nNode,1);
if param.flag2Dconstraints == 1
    rigid_mode_idx = [1 2 3];
else
    rigid_mode_idx = [1 2 3 4 5 6];
end
nMode = endMode-startMode+1-numel(rigid_mode_idx);
eVal(rigid_mode_idx,:) = [];
eVec(:,:,rigid_mode_idx) = [];
for jjj=1:nMode
    MSF(:,1) = MSF(:,1) + param.KbT/eVal(jjj,1)*sum(eVec(:,1:3,jjj).^2,2);
end
nodalRMSF = sqrt(MSF);

end