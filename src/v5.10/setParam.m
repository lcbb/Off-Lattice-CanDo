function param = setParam(path, angleHJ,rHelix,dampingMass,nickSF,alignBendHJ,alignTwistHJ, jobInfo)

% setting parameters for FE model

param.meanAngleHJ = deg2rad(angleHJ); % in radian
param.distBP = jobInfo.axialRise * 10; % the distance between base pairs in "ang"
param.angBP = 360/jobInfo.angBP; % 10.5 bp/turn, the twisting angle between base pairs in "deg"
param.rHelix = rHelix; % [Angstrom]
param.nu = 0.0;
param.A = pi*(param.rHelix)^2;
param.E = jobInfo.axialStiffness/param.A; % the effective Young's modulus in pN/ang^2
param.G = param.E/(2*(1+param.nu));
param.I = jobInfo.bendingStiffness*1e2/param.E; % the effective bending moment of inertia  in ang^4
param.J = jobInfo.torsionalStiffness*1e2/param.G; % the effective torsional moment of inertia in ang^4
param.rho = (2*326*0.1660538782e-4)/(param.A*param.distBP); % density in pN*tsec^2/ang^4
param.dampingMass = dampingMass; param.dampingStiff = 0.0; % uniform mass damping
param.nickSF = nickSF;
param.bulgeSF = nickSF * 1e-2;
param.alignStretchBP = 1e3*param.E*param.A/param.distBP;
param.alignBendBP = 1e3*param.E*param.I/param.distBP;
param.alignTwistBP = 1e3*param.G*param.J/param.distBP;
param.alignStretchHJ = 1e3*param.E*param.A/param.distBP;
param.alignBendHJ = alignBendHJ*param.E*param.I/param.distBP;
param.alignTwistHJ = alignTwistHJ*param.G*param.J/param.distBP;
param.alignBendHJ_ss = param.alignBendHJ;       % single crossover
param.alignTwistHJ_ss = param.alignTwistHJ * 1e-1;
param.alignBendHJ_bb = param.alignBendHJ * 2e-2;     % backbone
param.alignTwistHJ_bb = param.alignTwistHJ * 2e-1;

%%%%%%%%%%
param.ssDNA_axialStiffness = 800;
param.ssDNA_KuhnLength     = 15;
param.ssDNA_contourLength  = 5.6;
%%%%%%%%%%

param.model    = jobInfo.model;
param.atomic   = jobInfo.atomic;
param.nmaMovie = jobInfo.NMA;

param.krctypeHJ = 'CONSTANT';
param.krctypeBP = 'CONSTANT';
param.KbT = 1.3806504e-23*1e22*298; % T=298K, KbT in pN*ang
param.startMode = 1;
param.endMode = 206;
param.timestep = [ ...
                  1000  0.1;...
                  100  1.0; ...
                  50  10.0; ...
                  10  100.0; ...
                  10  1000.0; ...
                  10  10000.0; ...
                  10  100000.0; ...
                  10  1000000.0; ...
                  10  10000000.0; ...
                  10  100000000.0; ...
                  10  1000000000.0; ...
                 ];
param.cumulTimeStep = cumsum(param.timestep(:,1));
param.adinaSaveTimeStep = [ ...
                           1    10  1100  10; ...
                           2  1105  1150   5; ...
                           3  1152  1160   2; ...
                           4  1161  1230   1; ...
                          ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THESE PARAMETERS CONTROL AUXILIARY PROGRAMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param.auiEXE = path.ADINA_AUI;
param.auiOPTION = '-b -adina -m 8gb';
param.adinaEXE = path.ADINA;
param.adinaOPTION = '-b -s -mm 8gb -t 2';
param.chimearMatrixFN = 'chimera-matrix-view.txt';
param.chimeraEXE = path.CHIMERA;
param.chimeraOPTION = '--silent --script';
% Set video parameters if you want to make movies (by default not activated)
%param.videomachEXE = '"C:\Program Files (x86)\VideoMach\videomach.exe"';
%param.framerate = 10;

end
