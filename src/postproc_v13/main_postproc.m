function [] = main_postproc(designPATH,workDIR,prefixFN, param)
% Input files, including a *.out file and *.por file.
% !!! CAUTION !!! Some new files will be created in the directories
% containing these *.out and *.por files. The existing files with the same
% name will be overlapped.
%
% postprocDIR = 'D:\Postproc_offlattice';
% prefixFN = '21_21_2x4_nick1';
% prefixFN = '20_21_2x4_nick1';
% prefixFN = '20_21_32x4_nick1';
% prefixFN = '20_21_32x4_nick0.01';
% prefixFN = '22_21_2x4_nick1';
% prefixFN = '22_21_32x4_nick1';
% prefixFN = '22_21_32x4_nick0.01';
% prefixFN = '41_42_2x4_nick1';
% prefixFN = '43_42_2x4_nick1';
% prefixFN = '41_42_32x4_nick1';
% prefixFN = '41_42_32x4_nick0.01';
% prefixFN = '43_42_32x4_nick1';
% prefixFN = '43_42_32x4_nick0.01';

% postprocDIR = 'D:\Postproc_offlattice_v2';
% prefixFN = '20_21_2x4_nick0.01';
% prefixFN = 'DX_v3_nick0.01';
% prefixFN = 'DX_v3_tiamat_nick0.01';


%workDIR = fullfile(postprocDIR,prefixFN);
outPATH = fullfile(workDIR, strcat(prefixFN,'.out'));
porPATH = fullfile(workDIR, strcat(prefixFN,'.por'));
matPATH = fullfile(workDIR, strcat(prefixFN,'_FEModel.mat'));
msgPATH = fullfile(workDIR, strcat(prefixFN,'.msg'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lastStep = 1230;
lastTime = 0.1111111170000E+11;

% Time steps
param.timeStep = [10:10:1100, 1105:5:1150, 1152:2:1160 1161:1:1230]';

% Alignment element
param.alignHJE_ds.KTC1 = 323529.411765;
param.alignHJE_ds.KTC2 = 323529.411765;
param.alignHJE_ds.KTC3 = 323529.411765;
param.alignHJE_ds.KRC1 = 6764.705882;
param.alignHJE_ds.KRC2 = 6764.705882;
param.alignHJE_ds.KRC3 = 676.470588;

param.alignHJE_ss.KTC1 = 323529.411765;
param.alignHJE_ss.KTC2 = 323529.411765;
param.alignHJE_ss.KTC3 = 323529.411765;
param.alignHJE_ss.KRC1 = 0.676471;
param.alignHJE_ss.KRC2 = 0.676471;
param.alignHJE_ss.KRC3 = 0.067647;

param.alignDSDNA.KTC1 = 323529.411765;
param.alignDSDNA.KTC2 = 323529.411765;
param.alignDSDNA.KTC3 = 323529.411765;
param.alignDSDNA.KRC1 = 6764705.882353;
param.alignDSDNA.KRC2 = 6764705.882353;
param.alignDSDNA.KRC3 = 13529411.764706;

param.alignBBE = 323529.411765 * 1e-5;

% Control
% param.nmaMovie = false;
param.loadstepMovie = false;
% param.atomic = true;
param.atomicMovie = false;

% Visualization
% param.RMSF_range = [0,10];          % Unit: Angstrom
param.RMSF_range = [0,100];          % Unit: Angstrom
param.strE_range = [0,0.1];         % Unit: kBT
param.twist_range = [-1,1];         % Unit: degree
param.rHelix = 2.25*10/2;           % Unit: Angstrom
param.framerate = 10;
param.skiploadstep = 1;

% NMA
param.T_NMA = 298;
param.KbT = 1.3806504e-23*1e22 * param.T_NMA; % T=298K, KbT in pN*ang
% param.nFrame = 36;
param.nFrame = 72;
param.nPeriod = 1;
param.startMode = 1;
param.endMode = 206;

param.flag2Dconstraints = 0; % 1: fix out-of-plane displacement (x-displ)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post-processing finite element calculation results as 
% ADINA *.out and *.por files
%     1. Get the solution shape
%     2. Get the property of each node/element
%     3. Visualization
% Visualize strain energies from ADINA *.out and *.por files
% Render the solution shape as BILD and PNG files
% Five new files will be created in the working directory:
%     - *.bild
%     - *_chimeraScr.py
%     - *_view<1,2,3>.png
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Step 0. Check convergence
status = checkConvergence(msgPATH);
if(status ~= 0)
    warning('The current job did not converge.');
    return;
end
toc;

% Step 1. Get the displacement of each node from the *.por file
% Two new files will be created in the working directory:
%     - *_getDispl.plo
%     - *_displ.txt
% tic;
fprintf('Get displacements...');
displ = por2displ(porPATH,matPATH,param);
save(fullfile(workDIR,strcat(prefixFN,'_Displ.mat')), 'displ');
fprintf('Done.\n');
toc;

% Step 2. Get the coordinates of each node in the solution shape
fprintf('Get the solution shape...');
[dNode,pNode,pAxes] = getSolutionShape(matPATH,displ);
fprintf('Done.\n');
toc;

% Step 3. Strain energy
% Step 3.1 Get the strain energy for each beam element
fprintf('Get strain energy...');
strE = getStrainEnergy(outPATH,matPATH,lastStep,param);
save(fullfile(workDIR,strcat(prefixFN,'_StrE.mat')), 'strE');
% Step 3.2 Visualize strain energy for each beam element
%solutionShape2BildPng(matPATH,pNode,'_deformedShape_strE',[],strE,param);
solutionShape2BildPng(matPATH,pNode,'_deformedShape_strE',[],strE,param.strE_range,[],'jet',param,[]);
fprintf('Done.\n');
toc;

% Step 4. Normal mode analysis
if(param.nmaMovie)
    % Step 4.1 Get the RMSF value, eigenvalue (amplitude), and eigenvector
    %          (shape) for each node
    fprintf('Run NMA...');
    [eVal,eVec,nodalRMSF] = runNMA(porPATH,matPATH,param);
    save(fullfile(workDIR,strcat(prefixFN,'_NMA.mat')), 'eVal','eVec','nodalRMSF');
    % Step 4.2 Visualize RMSF for each node
    % fprintf('Rendering the solution shape...');
    solutionShape2BildPng(matPATH,pNode,'_deformedShape_RMSF',nodalRMSF,[],param.RMSF_range,[],'red',param,[]);
    % Step 4.3 Render the RMSF movie
    generateRMSFmovie_v2(workDIR,prefixFN,param,matPATH,pNode,pAxes,eVal,eVec,nodalRMSF);
    fprintf('Done.\n');
    toc;
end

% Step 5. Local geometry
% Step 5.1 Get the triad (using 3DNA convention) for each node
fprintf('Get local geometry...');
triad = getTriad(matPATH,porPATH,lastTime);
save(fullfile(workDIR,strcat(prefixFN,'_DOF.mat')), 'dNode','pNode','triad');
% Step 5.2 Get the 6 DOFs (using 3DNA convention) for each beam element
beamDOF = getBeamDOF(matPATH,dNode,triad);
save(fullfile(workDIR,strcat(prefixFN,'_3DNA.mat')), 'beamDOF');
% Step 5.3 Visualize overtwisting/undertwisting for each beam element
solutionShape2BildPng(matPATH,pNode,'_deformedShape_twist',[],beamDOF.twist-360/10.5,param.twist_range,[],'orange_blue',param,[]);
fprintf('Done.\n');
toc;

% Step 6. Calculate the total mechanical free energy
fprintf('Get total mechanical free energy at the last time step...');
[totStrE,beamStrE,trussStrE,alignStrE_HJE_ds,alignStrE_HJE_ss,alignStrE_DSDNA] = getTotalStrE_finalStep(matPATH,outPATH,param);
save(fullfile(workDIR,strcat(prefixFN,'_totStrE.mat')), 'totStrE','beamStrE','trussStrE','alignStrE_HJE_ds','alignStrE_HJE_ss','alignStrE_DSDNA');
fprintf('Done.\n');
fprintf('Total mechanical free energy: %.2f kBT\n', totStrE);
toc;


% Step 7. Analyze the load-step
if(param.loadstepMovie)
    % Step 7.1 Render the load-step movie
    fprintf('Rendering the load-step movie...');
    generateLoadStepmovie(workDIR,prefixFN,param,matPATH,outPATH);
    % Step 7.2 Plot the time history of mean nodal displacement magnitude and
    % total strain energy
    plotConvergence(workDIR,prefixFN,param);
    fprintf('Done.\n');
    toc;
end

% Step 8. Generate the full-atomistic model
if(param.atomic)
    main_generateAtomic(designPATH,workDIR,prefixFN,param);
end
toc
