function [] = main_offlattice(designs,model,atomic)

addpath src\postproc_v13
addpath src\DNA_topology_v17
addpath src\local_geometry_v2
addpath src\topology2pdb_ver4_bulge
addpath src\v5.10

% A. Input/output files
inputDIR  = 'input';
outputDIR = 'output';

% B. Setting Program Paths
path.ADINA_AUI = '"C:\ADINA93\x64\AUI.exe"';
path.ADINA     = '"C:\ADINA93\x64\adina.exe"';
path.CHIMERA   = '"C:\Program Files\Chimera 1.10.2\bin\chimera.exe"';

% C. Set Adjustable CanDo Parameters Here

% angleHJ [degree] : equilibrium junction twist angle (default = 60)
angleHJ = 60;

% rHelix [angstrom] : half of the interhelical distance (default = 0.5*18.5)
rHelix = 0.5*18.5;

% kTwist [pN nm rad^-1] : junction twist angle stiffness (default = 135.3)
kTwist = 135.3;

% nickSF [dimensionless] : nick stiffness factor (default = 1.0)
nickSF = 1.0e0;

% Other parameters which can be modified here are based on DNA duplex geometry
jobInfo.axialRise          = 0.34; % [nm]
jobInfo.angBP              = 10.5; % [bp / turn]
jobInfo.axialStiffness     = 1100; % [pN]
jobInfo.bendingStiffness   = 230;  % [pN nm2]
jobInfo.torsionalStiffness = 260;  % [pN nm2]
jobInfo.helixDiamater      = 2.25; % [nm]

% Parameter controlling the resolution of the model ('coarse' or 'fine' 
%(default)) 
jobInfo.model = model;

% Parameter controlling whether an atomic model is produced ('true' or 
%'false' (default))
jobInfo.atomic = atomic;

% Parameter controlling generation of normal mode analysis ('true' or 
%'false' (default))
jobInfo.NMA    = false;

% Other parameters - Be careful!!
dampingMass  = 1.0e0;    % default value
flagSSDNA    = 1;
alignBendHJ  = 1.0;
alignTwistHJ = 0.05;
tic;

% Main. Run FE simulations and Post-Processing
for iDesign=1:numel(designs)
                
    % Get current filenames & parameters
    bodyFN = designs{iDesign};
                
    % Create a working directory
    tarDIR = fullfile(outputDIR, bodyFN);
    if(~exist(tarDIR,'dir'))
        mkdir(tarDIR);
    end
                
    %% Setting CanDo parameters
    param = setParam(path,angleHJ,rHelix,dampingMass,nickSF,alignBendHJ,alignTwistHJ,jobInfo);

	%% Initialize type of topology file
	topTYPE = '';
	
    %% Parsing the input .dat, .cndo, and .txt files
    if ~exist(fullfile(inputDIR,strcat(bodyFN,'.dat'))) && ...,
            ~exist(fullfile(inputDIR,strcat(bodyFN,'.cndo')))
        topPATH = [];
        fprintf('Error: No topology file in input directory!\n');
    elseif exist(fullfile(inputDIR,strcat(bodyFN,'.dat')))
        topPATH = fullfile(inputDIR,strcat(bodyFN,'.dat'));
        fprintf('We have a Tiamat topology file. Setting the path...\n');
		topTYPE = 'tiamat';
	elseif exist(fullfile(inputDIR,strcat(bodyFN,'.cndo')))
        topPATH = fullfile(inputDIR,strcat(bodyFN,'.cndo'));
        fprintf('We have a CNDO topology file. Setting the path...\n');
		topTYPE = 'cndo';
    end

    if ~exist(fullfile(inputDIR,strcat(bodyFN,'.txt')))
        seqPATH = [];
        fprintf('Warning: No sequence file in input directory.\n');
    else
        seqPATH = fullfile(inputDIR,strcat(bodyFN,'.txt'));
        fprintf('We have a sequence file. Setting the path...\n');
    end

	if strcmp(topTYPE, 'tiamat')
	    fprintf('Processing %s ...', topPATH);
		main_Design2FE(topPATH, 'tiamat', '', outputDIR, seqPATH);
		fprintf('Done.\n');
		fprintf(1,'\t----- DNA Topology completed -----\n\n');
		toc;
	elseif strcmp(topTYPE, 'cndo')
	    fprintf('Processing %s ...', topPATH);
		main_Design2FE(topPATH, 'cndo', '', outputDIR, seqPATH);
		fprintf('Done.\n');
		fprintf(1,'\t----- DNA Topology completed -----\n\n');
		toc;
	else
		fprintf('No input design. Ending this iteration...\n');
		continue
	end
	
    % Locating structural motif information
    matHJE = fullfile(outputDIR,strcat(bodyFN,'.mat'));
                
    % Run off-lattice FE simulation
    fprintf('\n\n%s\n', repmat('=',1,80))
    fprintf('\tjob name=%s, nick stiffness=%.4f\n', bodyFN,nickSF);
    fprintf('%s\n\n', repmat('=',1,80))
    runFE_offlattice(matHJE,tarDIR,bodyFN,param,flagSSDNA);
                
    % Post-processing
    designPATH = fullfile(outputDIR,strcat(bodyFN,'.mat'));
    main_postproc(designPATH,tarDIR,bodyFN,param);

end
