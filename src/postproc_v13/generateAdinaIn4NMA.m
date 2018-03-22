function [] = generateAdinaIn4NMA(tarDIR,initFNbodyNLSA,initFNbodyNMA,param,FEModel)

% initFNbodyNLSA = strcat(bodyFN,param.analysisExtensionNonlinear);
% initFNbodyNMA = strcat(bodyFN,param.analysisExtensionNMA);

resFN_NLSA = fullfile(tarDIR,strcat(initFNbodyNLSA,'.res'));
idbFN_NLSA = fullfile(tarDIR,strcat(initFNbodyNLSA,'.idb'));

resFN_NMA = fullfile(tarDIR,strcat(initFNbodyNMA,'.res'));
idbFN_NMA = fullfile(tarDIR,strcat(initFNbodyNMA,'.idb'));
ilgFN_NMA = fullfile(tarDIR,strcat(initFNbodyNMA,'.ilg'));
datFN_NMA = fullfile(tarDIR,strcat(initFNbodyNMA,'.dat'));
inFN_NMA = fullfile(tarDIR,strcat(initFNbodyNMA,'.in'));

node = FEModel.node(:,1:3);
nNode = size(node,1);
nMode = min(param.endMode,nNode*6);

% copy .res file
if exist(resFN_NLSA,'file')~=2
    error('Restart file (%s) does not exist in %s',strcat(initFNbodyNLSA,'.res'),tarDIR);
else
    copyfile(resFN_NLSA,resFN_NMA,'f');
end


% generate .in file for NMA
fid = fopen(inFN_NMA,'w');
fprintf(fid,'FEPROGRAM ADINA\n');
fprintf(fid,'FILEECHO OPTION=FILE FILE=,\n');
fprintf(fid,'''');
fprintf(fid,'%s',ilgFN_NMA);
fprintf(fid,'''\n');
fprintf(fid,'FILELOG OPTION=FILE FILE=,\n');
fprintf(fid,'''');
fprintf(fid,'%s',ilgFN_NMA);
fprintf(fid,'''\n');
fprintf(fid,'DATABASE OPEN FILE=,\n');
fprintf(fid,'''');
fprintf(fid,'%s',idbFN_NLSA);
fprintf(fid,''' PERMFILE=,\n');
fprintf(fid,'''');
fprintf(fid,'%s',idbFN_NLSA);
fprintf(fid,''' PROMPT=NO\n');
fprintf(fid,'DELETE APPLY TYPE=DISPLACEMENTS SUBSTRUC=0 REUSE=1 THERMOST=0\n');
if param.flag2Dconstraints == 1
    load(fullfile(tarDIR,strcat(bodyFN,'.mat')),'caDNAnoData');
    node = caDNAnoData.node;
    fprintf(fid,'BOUNDARIES SUBSTRUC=0\n');
    fprintf(fid,'@CLEAR\n');
    for i=1:size(node,1)
        fprintf(fid,'%d ''FIXED''  ''FREE''  ''FREE''  ''FREE''  ''FREE''  ''FREE''  ''FREE''  ''FREE'',\n',i);
        fprintf(fid,'      ''FREE''  ''FREE''  ''FREE''  ''FREE''\n');
    end
    fprintf(fid,'9999999 ''FIXED''  ''FIXED''  ''FIXED''  ''FIXED''  ''FIXED''  ''FIXED''  ''FREE''  ''FREE'',\n');
    fprintf(fid,'      ''FREE''  ''FREE''  ''FREE''  ''FREE''\n');
    fprintf(fid,'@\n');
else
    fprintf(fid,'BOUNDARIES SUBSTRUC=0\n');
    fprintf(fid,'@CLEAR\n');
    fprintf(fid,'9999999 ''FIXED''  ''FIXED''  ''FIXED''  ''FIXED''  ''FIXED''  ''FIXED''  ''FREE''  ''FREE'',\n');
    fprintf(fid,'      ''FREE''  ''FREE''  ''FREE''  ''FREE''\n');
    fprintf(fid,'@\n');
end


fprintf(fid,'MASTER ANALYSIS=FREQUENCIES MODEX=RESTART TSTART=0.00000000000000,\n');
fprintf(fid,'     IDOF=0 OVALIZAT=NONE FLUIDPOT=AUTOMATIC CYCLICPA=1 IPOSIT=CONTINUE,\n');
fprintf(fid,'     REACTION=YES INITIALS=NO FSINTERA=NO IRINT=DEFAULT CMASS=NO,\n');
fprintf(fid,'     SHELLNDO=AUTOMATIC AUTOMATI=OFF SOLVER=SPARSE,\n');
fprintf(fid,'     CONTACT-=CONSTRAINT-FUNCTION TRELEASE=0.00000000000000,\n');
fprintf(fid,'     RESTART-=NO FRACTURE=NO LOAD-CAS=NO LOAD-PEN=NO SINGULAR=YES,\n');
fprintf(fid,'     STIFFNES=0.000100000000000000 MAP-OUTP=NONE MAP-FORM=NO,\n');
fprintf(fid,'     NODAL-DE='''' POROUS-C=NO ADAPTIVE=0 ZOOM-LAB=1 AXIS-CYC=0,\n');
fprintf(fid,'     PERIODIC=NO VECTOR-S=GEOMETRY EPSI-FIR=NO STABILIZ=NO,\n');
fprintf(fid,'     STABFACT=1.00000000000000E-10 RESULTS=PORTHOLE FEFCORR=NO,\n');
fprintf(fid,'     BOLTSTEP=1 EXTEND-S=YES CONVERT-=NO DEGEN=YES TMC-MODE=NO,\n');
fprintf(fid,'     ENSIGHT-=NO IRSTEPS=1 INITIALT=NO\n');
fprintf(fid,'FREQUENCIES METHOD=SUBSPACE-ITERATION NEIGEN=%d NMODE=%d IPRINT=NO,\n',nMode,nMode);
fprintf(fid,'     RIGID-BO=YES RSHIFT=-0.500000000000000,\n');
fprintf(fid,'     CUTOFF=1.00000000000000E+30 NITEMM=10000 NVECTOR=%d STURM-CH=NO,\n',2*nMode);
fprintf(fid,'     ACCELERA=NO TOLERANC=DEFAULT STARTTYP=LANCZOS NSTVECTO=0,\n');
fprintf(fid,'     INTERVAL=NO FMIN=0.00000000000000 FMAX=DEFAULT MODALSTR=NO,\n');
fprintf(fid,'     STATIC=NO NSHIFT=AUTO NSHIFT-B=50\n');
% fprintf(fid,'pproc nproc=4 minel=100\n');
fprintf(fid,'DATABASE SAVE PERMFILE=,\n');
fprintf(fid,'''');
fprintf(fid,'%s',idbFN_NMA);
fprintf(fid,''',\n');
fprintf(fid,'     PROMPT=NO\n');
fprintf(fid,'ADINA OPTIMIZE=SOLVER FILE=,\n');
fprintf(fid,'''');
fprintf(fid,'%s',datFN_NMA);
fprintf(fid,''',\n');
fprintf(fid,'     FIXBOUND=YES MIDNODE=NO OVERWRIT=YES\n');
fprintf(fid,'EXIT SAVE=NO IMMEDIATE=YES\n');
fclose(fid);
