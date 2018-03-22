function [] = generateRMSFmovie_v2(workDIR,prefixFN,param,matPATH,pNode,pAxes,eVal,eVec,nodalRMSF)

movieDIR = 'movieRMSF';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants and parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
startMode = param.startMode;      % main_CanDoWebServerVer2_proc1.m, Ln 369
endMode = param.endMode;      % main_CanDoWebServerVer2_proc1.m, Ln 370
rigid_mode_idx = [1 2 3 4 5 6];     % renderNMABildFormat.m, Ln 44

nNode = size(pNode,1);
endMode = min(endMode,nNode*6);

nMode = endMode-startMode+1-numel(rigid_mode_idx);
% eVal(rigid_mode_idx,:) = [];
% eVec(:,:,rigid_mode_idx) = [];

kBT = param.KbT;
nFrame = param.nFrame;
nPeriod = param.nPeriod;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the displacement in 6 DOFs in each frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi = zeros(nNode,6,nFrame);
tmax = nPeriod*(2*pi/sqrt(eVal(1,1)));
dt = tmax/nFrame;
t = dt:dt:tmax;

p0 = 2*pi*rand(nMode,1);
for i = 1:nFrame
    for j = 1:nMode
        phi(:,:,i) = phi(:,:,i) + sqrt(2*kBT/eVal(j))*eVec(:,:,j)*cos(sqrt(eVal(j,1))*t(i)+p0(j));
    end
end

save(fullfile(workDIR,strcat(prefixFN,'_phi.mat')), 'phi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pNodeSnap = zeros(nNode,3,nFrame);
for i = 1:nFrame
    pNodeSnap(:,:,i) = pNode + phi(:,1:3,i)*pAxes;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rendering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data array ---> BILD & PNG
solutionShape2BildPng(matPATH,pNodeSnap,'_RMSF',repmat(nodalRMSF,1,nFrame),[],param.RMSF_range,[],'red',param,movieDIR);

% PNG ---> AVI
% png2avi(fullfile(workDIR,movieDIR,strcat(prefixFN,'_RMSF')), size(pNodeSnap,3), param.framerate, param);

end