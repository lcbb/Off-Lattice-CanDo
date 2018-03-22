function [] = runFE_offlattice(matHJE,tarDIR,bodyFN,param,flagSSDNA)

a_groove = deg2rad(30);

% Loading structural motifs
load(matHJE, 'HJE','connectivity','nick','bulge','connSSE','connBBE','n_bulge');


% defining Holliday junction elements ------------------------------
elemHJE = struct([]);
elemDSE = struct([]);
nElemHJE = 0;
nElemHJE_ds = 0;
nElemHJE_ss = 0;
nElemDSE = 0;
for i=1:numel(HJE)
    if HJE(i).nWay == 2 % dsDNA segment
        nElemDSE = nElemDSE + 1;
        % [nm] in Tiamat --> [angstrom] in ADINA
        elemDSE(nElemDSE).cPos = 10*HJE(i).c'; % starting position in the e1 axis
        elemDSE(nElemDSE).cAxes = [HJE(i).e1, HJE(i).e2, HJE(i).e3]';
        elemDSE(nElemDSE).nBP = HJE(i).armLength;
        elemDSE(nElemDSE).idHJE = i;
        if isfield(HJE,{'isScaf'})
            elemDSE(nElemDSE).isScaf = HJE(i).isScaf;
        else
            elemDSE(nElemDSE).isScaf = false;
        end
    elseif HJE(i).nWay == 4 % Holliday junction element
        nElemHJE = nElemHJE + 1;
        if(strcmp(HJE(i).type, 'ds'))
            nElemHJE_ds = nElemHJE_ds + 1;
        elseif(strcmp(HJE(i).type, 'ss'))
            nElemHJE_ss = nElemHJE_ss + 1;
        else
            error('Exception.');
        end
        % [nm] in Tiamat --> [angstrom] in ADINA
        elemHJE(nElemHJE).cPos = 10*HJE(i).c';
        elemHJE(nElemHJE).cAxes = [HJE(i).e1, HJE(i).e2, HJE(i).e3]';
        elemHJE(nElemHJE).nBP_TLH = HJE(i).armLength(4); % number of basepairs in the top-left helix from the center
        elemHJE(nElemHJE).nBP_TRH = HJE(i).armLength(3); % number of basepairs in the top-right helix from the center
        elemHJE(nElemHJE).nBP_BLH = HJE(i).armLength(1); % number of basepairs in the bottom-left helix from the center
        elemHJE(nElemHJE).nBP_BRH = HJE(i).armLength(2); % number of basepairs in the bottom-right helix from the center
        elemHJE(nElemHJE).nBP = sum(HJE(i).armLength);
        elemHJE(nElemHJE).idHJE = i;
        elemHJE(nElemHJE).type = HJE(i).type;
        if isfield(HJE,{'isScaf'})
            elemHJE(nElemHJE).isScaf = HJE(i).isScaf;
        else
            elemHJE(nElemHJE).isScaf = false;
        end
    else
        error('Currently %d-way junctions cannot be modeled!\n',HJE(i).nWay);
    end
    assert(nElemHJE == nElemHJE_ds + nElemHJE_ss);
end

if isempty(elemDSE)
    if isempty(elemHJE)
        error('Neither HJE nor DSE exists!');
    else
        nNode = sum([elemHJE.nBP]);
    end
else
    if isempty(elemHJE)
        nNode = sum([elemDSE.nBP]);
    else
        nNode = sum([elemHJE.nBP])+sum([elemDSE.nBP]); % number of nodes = number of basepairs
    end
end
nElemDSDNA = nNode - 2*nElemHJE - nElemDSE;
%-------------------------------------------------------------------


% defining connectivity between Holliday junction elements ---------
% connHJE = [ ... eidHJ1, nidHJ1,  eidHJ2, nidHJ2 ... ]
% connSSE = [ ... eidHJ1, nidHJ1,  eidHJ2, nidHJ2, nBase ... ]
% nidHJ = 1(tip of BLH), 2(tip of BRH), 3(tip of TRH), 4(tip of TLH)
if(~isempty(connectivity))
    connHJE = connectivity(connectivity(:,5)==0&connectivity(:,6)==0,[1:4 7]);
    connHJEscaf = connHJE(connHJE(:,5)==1,1:2);
%     if flagSSDNA
%         connSSE = [connectivity(connectivity(:,5)~=0&connectivity(:,5)~=Inf,[1:4 5]); ...
%             connectivity(connectivity(:,6)~=0&connectivity(:,6)~=Inf,[1:4 6])];
%     else
%         connSSE = [];
%     end
else
    connHJE = [];
    connHJEscaf = zeros(0,2);    % connHJEscaf = [];
%     connSSE = [];
end
%-------------------------------------------------------------------

for i=1:nElemHJE
    % node number goes from left to right in the bottom helix
    %   nidT1----->-----nidC1-nidC2----->-----nidT2
    % node number goes from right to left in the top helix
    %   nidT2-----<-----nidC2-nidC1-----<-----nidT1
    if i==1
        elemHJE(1).nidT1_BH = 1;
        elemHJE(1).nidC1_BH = elemHJE(1).nidT1_BH + elemHJE(1).nBP_BLH - 1;
        elemHJE(1).nidC2_BH = elemHJE(1).nidC1_BH + 1;
        elemHJE(1).nidT2_BH = elemHJE(1).nidC2_BH + elemHJE(1).nBP_BRH - 1;
        elemHJE(1).nidT1_TH = elemHJE(1).nidT2_BH + 1;
        elemHJE(1).nidC1_TH = elemHJE(1).nidT1_TH + elemHJE(1).nBP_TRH - 1;
        elemHJE(1).nidC2_TH = elemHJE(1).nidC1_TH + 1;
        elemHJE(1).nidT2_TH = elemHJE(1).nidC2_TH + elemHJE(1).nBP_TLH - 1;
    else
        elemHJE(i).nidT1_BH = elemHJE(i-1).nidT2_TH + 1;
        elemHJE(i).nidC1_BH = elemHJE(i).nidT1_BH + elemHJE(i).nBP_BLH - 1;
        elemHJE(i).nidC2_BH = elemHJE(i).nidC1_BH + 1;
        elemHJE(i).nidT2_BH = elemHJE(i).nidC2_BH + elemHJE(i).nBP_BRH - 1;
        elemHJE(i).nidT1_TH = elemHJE(i).nidT2_BH + 1;
        elemHJE(i).nidC1_TH = elemHJE(i).nidT1_TH + elemHJE(i).nBP_TRH - 1;
        elemHJE(i).nidC2_TH = elemHJE(i).nidC1_TH + 1;
        elemHJE(i).nidT2_TH = elemHJE(i).nidC2_TH + elemHJE(i).nBP_TLH - 1;
    end
end

for i=1:nElemDSE
    % node number goes from left to right in the bottom helix
    %   nidT1----->-----nidT2
    if i==1
        if(nElemHJE == 0)
            elemDSE(1).nidT1 = 1;
        else
            elemDSE(1).nidT1 = elemHJE(nElemHJE).nidT2_TH + 1;
        end
        elemDSE(1).nidT2 = elemDSE(1).nidT1 + elemDSE(1).nBP - 1;
    else
        elemDSE(i).nidT1 = elemDSE(i-1).nidT2 + 1;
        elemDSE(i).nidT2 = elemDSE(i).nidT1 + elemDSE(i).nBP - 1;
    end
end

% generating FE nodes & beams for dsDNA
node = zeros(nNode,3);
nodalTriad = zeros(nNode,3,2);
elemDSDNA = zeros(nElemDSDNA,2);
nElemTMP = 0;
for i=1:nElemHJE
    
    if exist('nick','var')
        if ~isempty(nick)
            nickPos = nick(nick(:,1)==elemHJE(i).idHJE & nick(:,3)==0, 2);
            nickPos_open = nick(nick(:,1)==elemHJE(i).idHJE & nick(:,3)>0, 2);
        else
            nickPos = [];
            nickPos_open = [];
        end
    else
        nickPos = [];
        nickPos_open = [];
    end
    
    if exist('bulge','var')
        if ~isempty(bulge)
            bulgePos = bulge(bulge(:,1)==elemHJE(i).idHJE,2);
        else
            bulgePos = [];
        end
    else
        bulgePos = [];
    end
    
    % bottom helix
    dirBH = elemHJE(i).cAxes(1,:); % director vector of the helical axis (x-axis of HJ element unit)
    dirBH = dirBH/norm(dirBH);
    elemHJE(i).dirBH1 = dirBH;
    elemHJE(i).dirBH2 = elemHJE(i).cAxes(2,:)/norm(elemHJE(i).cAxes(2,:));
    elemHJE(i).dirBH3 = cross(elemHJE(i).dirBH1,elemHJE(i).dirBH2); elemHJE(i).dirBH3 = elemHJE(i).dirBH3/norm(elemHJE(i).dirBH3);
    cPosBH = elemHJE(i).cPos - param.rHelix*elemHJE(i).cAxes(3,:); % center position of the bottom helix
    distBH = param.distBP*( (-elemHJE(i).nBP_BLH:(elemHJE(i).nBP_BRH-1)) + 0.5)';
    angleBH = param.angBP*( (-elemHJE(i).nBP_BLH:(elemHJE(i).nBP_BRH-1)) + 0.5)'; % [degree]
    node(elemHJE(i).nidT1_BH:elemHJE(i).nidT2_BH,:) = repmat(cPosBH,elemHJE(i).nBP_BLH+elemHJE(i).nBP_BRH,1) + distBH*dirBH;
    nodalTriad(elemHJE(i).nidT1_BH:elemHJE(i).nidT2_BH,:,1) = repmat(elemHJE(i).dirBH1,elemHJE(i).nBP_BLH+elemHJE(i).nBP_BRH,1);
    nodalTriad(elemHJE(i).nidT1_BH:elemHJE(i).nidT2_BH,:,2) = cos(deg2rad(angleBH))*elemHJE(i).dirBH2 + sin(deg2rad(angleBH))*elemHJE(i).dirBH3;
    % elemDSDNA = [nid1, nid2, beam_type]
    % beam_type: 0 (dsDNA), 1 (nick)
    elemHJE(i).eidT1_BH = nElemTMP+1;
    elemHJE(i).eidT2_BH = nElemTMP+elemHJE(i).nBP_BLH+elemHJE(i).nBP_BRH-1;
    elemDSDNA(elemHJE(i).eidT1_BH:elemHJE(i).eidT2_BH,1) = (elemHJE(i).nidT1_BH:elemHJE(i).nidT2_BH-1)';
    elemDSDNA(elemHJE(i).eidT1_BH:elemHJE(i).eidT2_BH,2) = (elemHJE(i).nidT1_BH+1:elemHJE(i).nidT2_BH)';
    elemDSDNA(elemHJE(i).eidT1_BH:elemHJE(i).eidT2_BH,3) = zeros(numel(elemHJE(i).nidT1_BH+1:elemHJE(i).nidT2_BH),1);
    if ~isempty(nickPos)
        nickPosBH = nickPos( nickPos < (elemHJE(i).nBP_BLH+elemHJE(i).nBP_BRH) );
        if ~isempty(nickPosBH)
            elemDSDNA(nickPosBH+elemHJE(i).eidT1_BH-1,3) = 1; % beam_type = nick
        end
    end
    if ~isempty(nickPos_open)
        nickPosBH_open = nickPos_open( nickPos_open < (elemHJE(i).nBP_BLH+elemHJE(i).nBP_BRH) );
        elemDSDNA(nickPosBH_open+elemHJE(i).eidT1_BH-1,3) = 3; % beam_type = none
    end
    if ~isempty(bulgePos)
        bulgePosBH = bulgePos( bulgePos < (elemHJE(i).nBP_BLH+elemHJE(i).nBP_BRH) );
        if ~isempty(bulgePosBH)
            elemDSDNA(bulgePosBH+elemHJE(i).eidT1_BH-1,3) = 2; % beam_type = bulge
        end
    end
    nElemTMP = nElemTMP + elemHJE(i).nBP_BLH + elemHJE(i).nBP_BRH - 1;
    
    % top helix
    dirTH = cos(param.meanAngleHJ)*elemHJE(i).cAxes(1,:) + sin(param.meanAngleHJ)*elemHJE(i).cAxes(2,:);
    dirTH = -dirTH/norm(dirTH); % right to left
    elemHJE(i).dirTH1 = dirTH;
    elemHJE(i).dirTH2 = -sin(param.meanAngleHJ)*elemHJE(i).cAxes(1,:) + cos(param.meanAngleHJ)*elemHJE(i).cAxes(2,:);
    elemHJE(i).dirTH2 = -elemHJE(i).dirTH2/norm(elemHJE(i).dirTH2);
    elemHJE(i).dirTH3 = cross(elemHJE(i).dirTH1,elemHJE(i).dirTH2); elemHJE(i).dirTH3 = elemHJE(i).dirTH3/norm(elemHJE(i).dirTH3);
    cPosTH = elemHJE(i).cPos + param.rHelix*elemHJE(i).cAxes(3,:); % center position of the top helix
    distTH = param.distBP*( (-elemHJE(i).nBP_TRH:(elemHJE(i).nBP_TLH-1)) + 0.5)';
    angleTH = param.angBP*( (-elemHJE(i).nBP_TRH:(elemHJE(i).nBP_TLH-1)) + 0.5)'; % [degree]
    node(elemHJE(i).nidT1_TH:elemHJE(i).nidT2_TH,:) = repmat(cPosTH,elemHJE(i).nBP_TLH+elemHJE(i).nBP_TRH,1) + distTH*dirTH;
    nodalTriad(elemHJE(i).nidT1_TH:elemHJE(i).nidT2_TH,:,1) = repmat(elemHJE(i).dirTH1,elemHJE(i).nBP_TLH+elemHJE(i).nBP_TRH,1);
    nodalTriad(elemHJE(i).nidT1_TH:elemHJE(i).nidT2_TH,:,2) = -cos(deg2rad(angleTH))*elemHJE(i).dirTH2 - sin(deg2rad(angleTH))*elemHJE(i).dirTH3;
    % elemDSDNA = [nid1, nid2, beam_type]
    % beam_type: 0 (dsDNA), 1 (nick), 2 (bulge)
    elemHJE(i).eidT1_TH = nElemTMP+1;
    elemHJE(i).eidT2_TH = nElemTMP+elemHJE(i).nBP_TLH+elemHJE(i).nBP_TRH-1;
    elemDSDNA(elemHJE(i).eidT1_TH:elemHJE(i).eidT2_TH,1) = (elemHJE(i).nidT1_TH:elemHJE(i).nidT2_TH-1)';
    elemDSDNA(elemHJE(i).eidT1_TH:elemHJE(i).eidT2_TH,2) = (elemHJE(i).nidT1_TH+1:elemHJE(i).nidT2_TH)';
    elemDSDNA(elemHJE(i).eidT1_TH:elemHJE(i).eidT2_TH,3) = zeros(numel(elemHJE(i).nidT1_TH+1:elemHJE(i).nidT2_TH),1);
    if ~isempty(nickPos)
        nickPosTH = nickPos( nickPos > (elemHJE(i).nBP_BLH+elemHJE(i).nBP_BRH) );
        if ~isempty(nickPosTH)
            %                         elemDSDNA(nickPosTH+elemHJE(i).eidT1_TH-1,3) = 1; % beam_type = nick
            elemDSDNA(nickPosTH+elemHJE(i).eidT1_BH-2,3) = 1; % beam_type = nick
        end
    end
    if ~isempty(nickPos_open)
        nickPosTH_open = nickPos_open( nickPos_open > (elemHJE(i).nBP_BLH+elemHJE(i).nBP_BRH) );
        %                         elemDSDNA(nickPosTH+elemHJE(i).eidT1_TH-1,3) = 1; % beam_type = nick
        elemDSDNA(nickPosTH_open+elemHJE(i).eidT1_BH-2,3) = 3; % beam_type = none
    end
    if ~isempty(bulgePos)
        bulgePosTH = bulgePos( bulgePos > (elemHJE(i).nBP_BLH+elemHJE(i).nBP_BRH) );
        if ~isempty(bulgePosTH)
            %                         elemDSDNA(nickPosTH+elemHJE(i).eidT1_TH-1,3) = 1; % beam_type = nick
            elemDSDNA(bulgePosTH+elemHJE(i).eidT1_BH-2,3) = 2; % beam_type = bulge
        end
    end
    nElemTMP = nElemTMP + elemHJE(i).nBP_TLH + elemHJE(i).nBP_TRH - 1;
end


for i=1:nElemDSE
    
    nickPos = nick(nick(:,1)==elemDSE(i).idHJE & nick(:,3)==0, 2);
    nickPos_open = nick(nick(:,1)==elemDSE(i).idHJE & nick(:,3)>0, 2);
    
    % bottom helix (no top helix, in fact)
    dirBH = elemDSE(i).cAxes(1,:); % director vector of the helical axis (x-axis of HJ element unit)
    dirBH = dirBH/norm(dirBH);
    elemDSE(i).dirBH1 = dirBH;
    elemDSE(i).dirBH2 = elemDSE(i).cAxes(2,:)/norm(elemDSE(i).cAxes(2,:));
    elemDSE(i).dirBH3 = cross(elemDSE(i).dirBH1,elemDSE(i).dirBH2); elemDSE(i).dirBH3 = elemDSE(i).dirBH3/norm(elemDSE(i).dirBH3);
    cPosBH = elemDSE(i).cPos; % starting basepair position
    distBH = param.distBP*( 0:(elemDSE(i).nBP-1) )';
    angleBH = param.angBP*( 0:(elemDSE(i).nBP-1) )'; % [degree]
    node(elemDSE(i).nidT1:elemDSE(i).nidT2,:) = repmat(cPosBH,elemDSE(i).nBP,1) + distBH*dirBH;
    nodalTriad(elemDSE(i).nidT1:elemDSE(i).nidT2,:,1) = repmat(elemDSE(i).dirBH1,elemDSE(i).nBP,1);
    nodalTriad(elemDSE(i).nidT1:elemDSE(i).nidT2,:,2) = cos(deg2rad(angleBH))*elemDSE(i).dirBH2 + sin(deg2rad(angleBH))*elemDSE(i).dirBH3;
    % elemDSDNA = [nid1, nid2, beam_type]
    % beam_type: 0 (dsDNA), 1 (nick)
    elemDSE(i).eidT1 = nElemTMP+1;
    elemDSE(i).eidT2 = nElemTMP+elemDSE(i).nBP-1;
    elemDSDNA(elemDSE(i).eidT1:elemDSE(i).eidT2,1) = (elemDSE(i).nidT1:elemDSE(i).nidT2-1)';
    elemDSDNA(elemDSE(i).eidT1:elemDSE(i).eidT2,2) = (elemDSE(i).nidT1+1:elemDSE(i).nidT2)';
    elemDSDNA(elemDSE(i).eidT1:elemDSE(i).eidT2,3) = zeros(numel(elemDSE(i).nidT1+1:elemDSE(i).nidT2),1);
    if ~isempty(nickPos)
        nickPosBH = nickPos( nickPos < (elemDSE(i).nBP) );
        if ~isempty(nickPosBH)
            elemDSDNA(nickPosBH+elemDSE(i).eidT1-1,3) = 1; % beam_type = nick
        end
    end
    if ~isempty(nickPos_open)
        nickPosBH_open = nickPos_open( nickPos_open < (elemDSE(i).nBP) );
        elemDSDNA(nickPosBH_open+elemDSE(i).eidT1-1,3) = 3; % beam_type = none
    end
    nElemTMP = nElemTMP + elemDSE(i).nBP - 1;
end


% writing ADINA input file ---------
nu = param.nu; E = param.E; rho = param.rho; II = param.I; JJ = param.J; AA = param.A;

fid = fopen(fullfile(tarDIR,strcat(bodyFN,'.in')),'w');
fprintf(fid,'DATABASE NEW SAVE=NO PROMPT=NO\n');
fprintf(fid,'FEPROGRAM ADINA\n');
fprintf(fid,'CONTROL PLOTUNIT=PERCENT VERBOSE=YES ERRORLIM=0 LOGLIMIT=0 UNDO=-1,\n');
fprintf(fid,'     PROMPTDE=NO AUTOREPA=NO DRAWMATT=YES DRAWTEXT=EXACT,\n');
fprintf(fid,'     DRAWLINE=EXACT DRAWFILL=EXACT AUTOMREB=NO ZONECOPY=NO,\n');
fprintf(fid,'     SWEEPCOI=YES SESSIONS=YES DYNAMICT=YES UPDATETH=YES AUTOREGE=NO,\n');
fprintf(fid,'     ERRORACT=CONTINUE FILEVERS=V88 INITFCHE=NO SIGDIGIT=6,\n');
fprintf(fid,'     AUTOZONE=YES PSFILEVE=V0\n');
fprintf(fid,'KINEMATICS DISPLACE=LARGE STRAINS=SMALL UL-FORMU=DEFAULT PRESSURE=NO,\n');
fprintf(fid,'     INCOMPAT=AUTOMATIC RIGIDLIN=NO BEAM-ALGORITHM=CURRENT\n');
fprintf(fid,'PORTHOLE SAVEDEFA=YES FILEUNIT=60 FORMATTE=YES INPUT-DA=1,\n');
fprintf(fid,'     DISPLACE=YES VELOCITI=YES ACCELERA=YES TEMPERAT=YES MAX-STEP=0,\n');
fprintf(fid,'     SHELLVEC=NO ELEM-RES=DEF-GRAD\n');
fprintf(fid,'****ADD:BATCH\n');
fprintf(fid,'FILEECHO OPTION=FILE FILE=''%s''\n',fullfile(tarDIR,strcat(bodyFN,'.ilg')));
fprintf(fid,'FILELOG OPTION=FILE FILE=''%s''\n',fullfile(tarDIR,strcat(bodyFN,'.ilg')));
fprintf(fid,'****END:BATCH\n');

% time step
fprintf(fid,'TIMESTEP NAME=DEFAULT\n');
fprintf(fid,'@CLEAR\n');
for i=1:size(param.timestep,1)
    fprintf(fid,'%d %f\n',param.timestep(i,1),param.timestep(i,2));
end
fprintf(fid,'@\n');

% coordinates
FEModel.node = node;
fprintf(fid,'COORDINATES NODE SYSTEM=0\n');
fprintf(fid,'@CLEAR\n');
for i=1:size(node,1)
    fprintf(fid,'%d  %f  %f  %f  0\n',i,FEModel.node(i,1),FEModel.node(i,2),FEModel.node(i,3));
end
fprintf(fid,'9999999  %f  %f  %f  0\n',FEModel.node(i,1)-1.013,FEModel.node(i,2)+1.0123,FEModel.node(i,3)-0.1231);
fprintf(fid,'@\n');

% boundary condition
[~,centerNodeID] = min(sum((node - repmat(mean(node,1),nNode,1)).^2,2));
fprintf(fid,'BOUNDARIES SUBSTRUC=0\n');
fprintf(fid,'@CLEAR\n');
fprintf(fid,'%d ''FIXED''  ''FIXED''  ''FIXED''  ''FIXED''  ''FIXED''  ''FIXED''  ''FREE''  ''FREE'',\n',centerNodeID);
fprintf(fid,'      ''FREE''  ''FREE''  ''FREE''  ''FREE''\n');
fprintf(fid,'9999999 ''FIXED''  ''FIXED''  ''FIXED''  ''FIXED''  ''FIXED''  ''FIXED''  ''FREE''  ''FREE'',\n');
fprintf(fid,'      ''FREE''  ''FREE''  ''FREE''  ''FREE''\n');
fprintf(fid,'@\n');

nElemDSDNA = size(elemDSDNA,1);
normalDSDNA = elemDSDNA(elemDSDNA(:,3)==0,1:2);
nickedDSDNA = elemDSDNA(elemDSDNA(:,3)==1,1:2);
bulgedDSDNA = elemDSDNA(elemDSDNA(:,3)==2,1:2);
noneDSDNA = elemDSDNA(elemDSDNA(:,3)==3,1:2);
nNormalDSDNA = size(normalDSDNA,1);
nNickedDSDNA = size(nickedDSDNA,1);
nBulgedDSDNA = size(bulgedDSDNA,1);
nNoneDSDNA = size(noneDSDNA,1);

% normal beams
if nNormalDSDNA ~=0
    FEModel.beamsNormal = zeros(nNormalDSDNA,3);
    FEModel.beamsNormal(:,1) = (1:nNormalDSDNA)';
    FEModel.beamsNormal(:,2:3) = normalDSDNA;
    fprintf(fid,'MATERIAL ELASTIC NAME=1 E=%f NU=%f,\n',E,nu);
    fprintf(fid,'     DENSITY=%f ALPHA=0.00000000000000 MDESCRIP=''NONE''\n',rho);
    fprintf(fid,'CROSS-SECTION PROPERTIES NAME=1,\n');
    fprintf(fid,'     RINERTIA=%f  SINERTIA=%f,\n',JJ,II);
    fprintf(fid,'     TINERTIA=%f  AREA=%f\n',II,AA);
    fprintf(fid,'EGROUP BEAM NAME=1 SUBTYPE=THREE-D DISPLACE=DEFAULT MATERIAL=1 RINT=5,\n');
    fprintf(fid,'     SINT=DEFAULT TINT=DEFAULT RESULTS=STRESSES INITIALS=NONE,\n');
    fprintf(fid,'     CMASS=DEFAULT RIGIDEND=NONE MOMENT-C=NO RIGIDITY=1,\n');
    fprintf(fid,'     MULTIPLY=1000000.00000000 RUPTURE=ADINA OPTION=NONE,\n');
    fprintf(fid,'     BOLT-TOL=0.00000000000000 DESCRIPT=''NONE'' SECTION=1,\n');
    fprintf(fid,'     PRINT=DEFAULT SAVE=DEFAULT TBIRTH=0.00000000000000,\n');
    fprintf(fid,'     TDEATH=0.00000000000000 SPOINT=4 BOLTFORC=0.00000000000000,\n');
    fprintf(fid,'     BOLTNCUR=0 TMC-MATE=1 BOLT-NUM=0 BOLT-LOA=0.00000000000000,\n');
    fprintf(fid,'     WARP=NO\n');
    fprintf(fid,'ENODES SUBSTRUC=0 GROUP=1 NNODES=32\n');
    fprintf(fid,'@CLEAR\n');
    for i=1:nNormalDSDNA
        fprintf(fid,'%d  9999999  %d  %d\n',FEModel.beamsNormal(i,1),FEModel.beamsNormal(i,2),FEModel.beamsNormal(i,3));
    end
    fprintf(fid,'@\n');
else
    FEModel.beamsNormal = [];
end

% nicked beams
if nNickedDSDNA ~=0
    FEModel.beamsNick = zeros(nNickedDSDNA,3);
    FEModel.beamsNick(:,1) = (1:nNickedDSDNA)'+nNormalDSDNA;
    FEModel.beamsNick(:,2:3) = nickedDSDNA;
    fprintf(fid,'MATERIAL ELASTIC NAME=2 E=%f NU=%f,\n',E,nu);
    fprintf(fid,'     DENSITY=%f ALPHA=0.00000000000000 MDESCRIP=''NONE''\n',rho);
    fprintf(fid,'CROSS-SECTION PROPERTIES NAME=2,\n');
    fprintf(fid,'     RINERTIA=%f  SINERTIA=%f,\n',JJ*param.nickSF,II*param.nickSF);
    fprintf(fid,'     TINERTIA=%f  AREA=%f\n',II*param.nickSF,AA);
    fprintf(fid,'EGROUP BEAM NAME=2 SUBTYPE=THREE-D DISPLACE=DEFAULT MATERIAL=2 RINT=5,\n');
    fprintf(fid,'     SINT=DEFAULT TINT=DEFAULT RESULTS=STRESSES INITIALS=NONE,\n');
    fprintf(fid,'     CMASS=DEFAULT RIGIDEND=NONE MOMENT-C=NO RIGIDITY=1,\n');
    fprintf(fid,'     MULTIPLY=1000000.00000000 RUPTURE=ADINA OPTION=NONE,\n');
    fprintf(fid,'     BOLT-TOL=0.00000000000000 DESCRIPT=''NONE'' SECTION=2,\n');
    fprintf(fid,'     PRINT=DEFAULT SAVE=DEFAULT TBIRTH=0.00000000000000,\n');
    fprintf(fid,'     TDEATH=0.00000000000000 SPOINT=4 BOLTFORC=0.00000000000000,\n');
    fprintf(fid,'     BOLTNCUR=0 TMC-MATE=1 BOLT-NUM=0 BOLT-LOA=0.00000000000000,\n');
    fprintf(fid,'     WARP=NO\n');
    fprintf(fid,'ENODES SUBSTRUC=0 GROUP=2 NNODES=32\n');
    fprintf(fid,'@CLEAR\n');
    for i=1:nNickedDSDNA
        fprintf(fid,'%d  9999999  %d  %d\n',FEModel.beamsNick(i,1),FEModel.beamsNick(i,2),FEModel.beamsNick(i,3));
    end
    fprintf(fid,'@\n');
else
    FEModel.beamsNick = [];
end

% bulged beams
if nBulgedDSDNA ~=0
    FEModel.beamsBulge = zeros(nBulgedDSDNA,3);
    FEModel.beamsBulge(:,1) = (1:nBulgedDSDNA)'+nNormalDSDNA+nNickedDSDNA;
    FEModel.beamsBulge(:,2:3) = bulgedDSDNA;
    fprintf(fid,'MATERIAL ELASTIC NAME=3 E=%f NU=%f,\n',E,nu);
    fprintf(fid,'     DENSITY=%f ALPHA=0.00000000000000 MDESCRIP=''NONE''\n',rho);
    fprintf(fid,'CROSS-SECTION PROPERTIES NAME=3,\n');
    fprintf(fid,'     RINERTIA=%f  SINERTIA=%f,\n',JJ*param.bulgeSF,II*param.bulgeSF);
    fprintf(fid,'     TINERTIA=%f  AREA=%f\n',II*param.bulgeSF,AA);
    fprintf(fid,'EGROUP BEAM NAME=3 SUBTYPE=THREE-D DISPLACE=DEFAULT MATERIAL=3 RINT=5,\n');
    fprintf(fid,'     SINT=DEFAULT TINT=DEFAULT RESULTS=STRESSES INITIALS=NONE,\n');
    fprintf(fid,'     CMASS=DEFAULT RIGIDEND=NONE MOMENT-C=NO RIGIDITY=1,\n');
    fprintf(fid,'     MULTIPLY=1000000.00000000 RUPTURE=ADINA OPTION=NONE,\n');
    fprintf(fid,'     BOLT-TOL=0.00000000000000 DESCRIPT=''NONE'' SECTION=3,\n');
    fprintf(fid,'     PRINT=DEFAULT SAVE=DEFAULT TBIRTH=0.00000000000000,\n');
    fprintf(fid,'     TDEATH=0.00000000000000 SPOINT=4 BOLTFORC=0.00000000000000,\n');
    fprintf(fid,'     BOLTNCUR=0 TMC-MATE=1 BOLT-NUM=0 BOLT-LOA=0.00000000000000,\n');
    fprintf(fid,'     WARP=NO\n');
    fprintf(fid,'ENODES SUBSTRUC=0 GROUP=3 NNODES=32\n');
    fprintf(fid,'@CLEAR\n');
    for i=1:nBulgedDSDNA
        fprintf(fid,'%d  9999999  %d  %d\n',FEModel.beamsBulge(i,1),FEModel.beamsBulge(i,2),FEModel.beamsBulge(i,3));
    end
    fprintf(fid,'@\n');
else
    FEModel.beamsBulge = [];
end

% none beams
if nNoneDSDNA ~=0
    FEModel.beamsNone = zeros(nNoneDSDNA,3);
    FEModel.beamsNone(:,1) = (1:nNoneDSDNA)'+nNormalDSDNA+nNickedDSDNA+nBulgedDSDNA;
    FEModel.beamsNone(:,2:3) = noneDSDNA;
    fprintf(fid,'MATERIAL ELASTIC NAME=4 E=%f NU=%f,\n',E,nu);
    fprintf(fid,'     DENSITY=%f ALPHA=0.00000000000000 MDESCRIP=''NONE''\n',rho);
    fprintf(fid,'CROSS-SECTION PROPERTIES NAME=4,\n');
    fprintf(fid,'     RINERTIA=%f  SINERTIA=%f,\n',JJ*1e-6,II*1e-6);
    fprintf(fid,'     TINERTIA=%f  AREA=%f\n',II*1e-6,AA*1e-6);
    fprintf(fid,'EGROUP BEAM NAME=4 SUBTYPE=THREE-D DISPLACE=DEFAULT MATERIAL=4 RINT=5,\n');
    fprintf(fid,'     SINT=DEFAULT TINT=DEFAULT RESULTS=STRESSES INITIALS=NONE,\n');
    fprintf(fid,'     CMASS=DEFAULT RIGIDEND=NONE MOMENT-C=NO RIGIDITY=1,\n');
    fprintf(fid,'     MULTIPLY=1000000.00000000 RUPTURE=ADINA OPTION=NONE,\n');
    fprintf(fid,'     BOLT-TOL=0.00000000000000 DESCRIPT=''NONE'' SECTION=4,\n');
    fprintf(fid,'     PRINT=DEFAULT SAVE=DEFAULT TBIRTH=0.00000000000000,\n');
    fprintf(fid,'     TDEATH=0.00000000000000 SPOINT=4 BOLTFORC=0.00000000000000,\n');
    fprintf(fid,'     BOLTNCUR=0 TMC-MATE=1 BOLT-NUM=0 BOLT-LOA=0.00000000000000,\n');
    fprintf(fid,'     WARP=NO\n');
    fprintf(fid,'ENODES SUBSTRUC=0 GROUP=4 NNODES=32\n');
    fprintf(fid,'@CLEAR\n');
    for i=1:nNoneDSDNA
        fprintf(fid,'%d  9999999  %d  %d\n',FEModel.beamsNone(i,1),FEModel.beamsNone(i,2),FEModel.beamsNone(i,3));
    end
    fprintf(fid,'@\n');
else
    FEModel.beamsNone = [];
end

% ssDNA connections
if ~isempty(connSSE)
    for i=1:size(connSSE,1)
        
        hjeid1 = connSSE(i,1);
        if HJE(hjeid1).nWay == 4
            idx1 = find([elemHJE.idHJE]==hjeid1,1);
            nid1 = elemHJE(idx1).nidT1_BH - 1 + connSSE(i,2);
%             switch(connSSE(i,2))
%                 case 1,
%                     nid1 = elemHJE(idx1).nidT1_BH;
%                 case 2,
%                     nid1 = elemHJE(idx1).nidT2_BH;
%                 case 3,
%                     nid1 = elemHJE(idx1).nidT1_TH;
%                 case 4,
%                     nid1 = elemHJE(idx1).nidT2_TH;
%             end
        elseif HJE(hjeid1).nWay == 2
            idx1 = find([elemDSE.idHJE]==hjeid1,1);
            nid1 = elemDSE(idx1).nidT1 - 1 + connSSE(i,2);
%             switch(connSSE(i,2))
%                 case 1,
%                     nid1 = elemDSE(idx1).nidT1;
%                 case 2,
%                     nid1 = elemDSE(idx1).nidT2;
%             end
        else
            error('Unknown %d-way junction!\n',HJE(hjeid1).nWay);
        end
        
        hjeid2 = connSSE(i,3);
        if HJE(hjeid2).nWay == 4
            idx2 = find([elemHJE.idHJE]==hjeid2,1);
            nid2 = elemHJE(idx2).nidT1_BH - 1 + connSSE(i,4);
%             switch(connSSE(i,4))
%                 case 1,
%                     nid2 = elemHJE(idx2).nidT1_BH;
%                 case 2
%                     nid2 = elemHJE(idx2).nidT2_BH;
%                 case 3,
%                     nid2 = elemHJE(idx2).nidT1_TH;
%                 case 4,
%                     nid2 = elemHJE(idx2).nidT2_TH;
%             end
        elseif HJE(hjeid2).nWay == 2
            idx2 = find([elemDSE.idHJE]==hjeid2,1);
            nid2 = elemDSE(idx2).nidT1 - 1 + connSSE(i,4);
%             switch(connSSE(i,4))
%                 case 1,
%                     nid2 = elemDSE(idx2).nidT1;
%                 case 2,
%                     nid2 = elemDSE(idx2).nidT2;
%             end
        else
            error('Unknown %d-way junction!\n',HJE(hjeid2).nWay);
        end
        
        nBase = connSSE(i,5);
        distNode = norm( node(nid1,1:3)-node(nid2,1:3) );
        %if(nBase >= 5)
        if(nBase >= 1)
            ss1 =-10:9:-1; ss2 = [1:10 15:5:50 60:10:500 1e20];
            ss1
            ss2
            %ee1 = (5*nBase)/distNode*(coth(ss1*15/param.KbT)-param.KbT./ss1/15).*(1+ss1/800) - 1.0;
            %ee2 = (5*nBase)/distNode*(coth(ss2*15/param.KbT)-param.KbT./ss2/15).*(1+ss2/800) - 1.0;
            ee1 = (param.ssDNA_contourLength*nBase)/distNode*(coth(ss1*param.ssDNA_KuhnLength/param.KbT)-param.KbT./ss1/param.ssDNA_KuhnLength).*(1+ss1/param.ssDNA_axialStiffness) - 1.0;
            ee2 = (param.ssDNA_contourLength*nBase)/distNode*(coth(ss2*param.ssDNA_KuhnLength/param.KbT)-param.KbT./ss2/param.ssDNA_KuhnLength).*(1+ss2/param.ssDNA_axialStiffness) - 1.0;
            ee = [-1000 ee1 -1.0 ee2];
            ee
            ss = [-1000 ss1  0.0 ss2];
            ss
        else
            ss = [-1000 -10:9:-1 0 1:10 15:5:50 60:10:500 1e20];
            ee = (3.4*nBase+3.4)/distNode.*(1+ss/800) - 1.0;
        end
        
        fprintf(fid,'MATERIAL NONLINEAR-ELASTIC NAME=%d DENSITY=%f,\n',i+4,rho);
        fprintf(fid,'     DCURVE=0 MDESCRIP=''NONE'' NU=0.0 MATRIX=TANGENT\n');
        fprintf(fid,'@CLEAR\n');
        for j=1:numel(ee)
            fprintf(fid,'%f  %f\n',ee(j),ss(j));
        end
        fprintf(fid,'@\n');
        fprintf(fid,'EGROUP TRUSS NAME=%d SUBTYPE=GENERAL DISPLACE=DEFAULT MATERIAL=%d,\n',i+4,i+4);
        fprintf(fid,'     INT=DEFAULT GAPS=NO INITIALS=ELEMENT CMASS=DEFAULT,\n');
        fprintf(fid,'     TIME-OFF=0.00000000000000 OPTION=NONE RB-LINE=1 DESCRIPT=''NONE'',\n');
        fprintf(fid,'     AREA=1.00000000000000 PRINT=DEFAULT SAVE=DEFAULT,\n');
        fprintf(fid,'     TBIRTH=0.00000000000000 TDEATH=0.00000000000000 TMC-MATE=1,\n');
        fprintf(fid,'     RUPTURE=ADINA GAPWIDTH=0.00000000000000\n');
        fprintf(fid,'ENODES SUBSTRUC=0 GROUP=%d NNODES=32\n',i+4);
        fprintf(fid,'@CLEAR\n');
        fprintf(fid,'%d  %d  %d  0  0\n',i+nNormalDSDNA+nNickedDSDNA+nBulgedDSDNA+nNoneDSDNA, ...
            nid1,nid2);
        fprintf(fid,'@\n');
        FEModel.entSprings(i,:) = [i+nNormalDSDNA+nNickedDSDNA+nBulgedDSDNA+nNoneDSDNA, nid1, nid2];
        
    end
    nEntSprings = size(connSSE,1);
    
else
    nEntSprings = 0;
    FEModel.entSprings = [];
end

% backbone connections
% format: each row = [node_id_1, triadset_id_1, node_id_2, triadset_id_2]
if ~isempty(connBBE)
    for i = 1:size(connBBE,1)
        % 5'-end of the backbone connection
        hjeid1 = connBBE(i,1);
        if HJE(hjeid1).nWay == 4
            idx1 = find([elemHJE.idHJE]==hjeid1,1);
            nid1 = elemHJE(idx1).nidT1_BH - 1 + connBBE(i,2);
            if(nid1 == elemHJE(idx1).nidT1_BH)
                triadset_id_1 = 16 * (idx1-1) + 7;
            elseif(nid1 == elemHJE(idx1).nidT2_BH)
                triadset_id_1 = 16 * (idx1-1) + 10;
            elseif(nid1 == elemHJE(idx1).nidT1_TH)
                triadset_id_1 = 16 * (idx1-1) + 13;
            elseif(nid1 == elemHJE(idx1).nidT2_TH)
                triadset_id_1 = 16 * (idx1-1) + 16;
            else
                error('Exception.');
            end
        elseif HJE(hjeid1).nWay == 2
            idx1 = find([elemDSE.idHJE]==hjeid1,1);
            nid1 = elemDSE(idx1).nidT1 - 1 + connBBE(i,2);
            if(nid1 == elemDSE(idx1).nidT1)
                triadset_id_1 = 16 * nElemHJE + 6 * (idx1-1) + 3;
            elseif(nid1 == elemDSE(idx1).nidT2)
                triadset_id_1 = 16 * nElemHJE + 6 * (idx1-1) + 6;
            end
        else
            error('Unknown %d-way junction!\n',HJE(hjeid1).nWay);
        end
        
        % 3'-end of the backbone connection
        hjeid2 = connBBE(i,3);
        if HJE(hjeid2).nWay == 4
            idx2 = find([elemHJE.idHJE]==hjeid2,1);
            nid2 = elemHJE(idx2).nidT1_BH - 1 + connBBE(i,4);
            if(nid2 == elemHJE(idx2).nidT1_BH)
                triadset_id_2 = 16 * (idx2-1) + 6;
            elseif(nid2 == elemHJE(idx2).nidT2_BH)
                triadset_id_2 = 16 * (idx2-1) + 9;
            elseif(nid2 == elemHJE(idx2).nidT1_TH)
                triadset_id_2 = 16 * (idx2-1) + 12;
            elseif(nid2 == elemHJE(idx2).nidT2_TH)
                triadset_id_2 = 16 * (idx2-1) + 15;
            else
                error('Exception.');
            end
        elseif HJE(hjeid2).nWay == 2
            idx2 = find([elemDSE.idHJE]==hjeid2,1);
            nid2 = elemDSE(idx2).nidT1 - 1 + connBBE(i,4);
            if(nid2 == elemDSE(idx2).nidT1)
                triadset_id_2 = 16 * nElemHJE + 6 * (idx2-1) + 2;
            elseif(nid2 == elemDSE(idx2).nidT2)
                triadset_id_2 = 16 * nElemHJE + 6 * (idx2-1) + 5;
            end
        else
            error('Unknown %d-way junction!\n',HJE(hjeid2).nWay);
        end
        
        FEModel.BBE(i,:) = [nid1, triadset_id_1, nid2, triadset_id_2];
    end
    nBBE = size(FEModel.BBE, 1);
else
    nBBE = 0;
    FEModel.BBE = [];
end

% alignment definition
fprintf(fid,'ALIGN-TRANSLATION  1  KTC1=%f  KTC2=%f  KTC3=%f\n',param.alignStretchHJ,param.alignStretchHJ,param.alignStretchHJ);
fprintf(fid,'@CLEAR\n');
fprintf(fid,'0   TBFACTOR\n');
fprintf(fid,'100  ALIGNED\n');
fprintf(fid,'@\n');
fprintf(fid,'ALIGN-ROTATION  1  KRC1=%f  KRC2=%f  KRC3=%f,\n',param.alignBendHJ,param.alignBendHJ,param.alignTwistHJ);
fprintf(fid,'  KRC1TYPE=%s  KRC2TYPE=%s  KRC3TYPE=%s\n',param.krctypeHJ,param.krctypeHJ,param.krctypeHJ);
fprintf(fid,'@CLEAR\n');
fprintf(fid,'0   TBFACTOR\n');
fprintf(fid,'100 ALIGNED\n');
fprintf(fid,'@\n');
fprintf(fid,'ALIGN-TRANSLATION  2  KTC1=%f  KTC2=%f  KTC3=%f\n',param.alignStretchHJ,param.alignStretchHJ,param.alignStretchHJ);
fprintf(fid,'@CLEAR\n');
fprintf(fid,'0   TBFACTOR\n');
fprintf(fid,'100  ALIGNED\n');
fprintf(fid,'@\n');
fprintf(fid,'ALIGN-ROTATION  2  KRC1=%f  KRC2=%f  KRC3=%f,\n',param.alignBendHJ_ss,param.alignBendHJ_ss,param.alignTwistHJ_ss);
fprintf(fid,'  KRC1TYPE=%s  KRC2TYPE=%s  KRC3TYPE=%s\n',param.krctypeHJ,param.krctypeHJ,param.krctypeHJ);
fprintf(fid,'@CLEAR\n');
fprintf(fid,'0   TBFACTOR\n');
fprintf(fid,'100 ALIGNED\n');
fprintf(fid,'@\n');
fprintf(fid,'ALIGN-TRANSLATION  3  KTC1=%f  KTC2=%f  KTC3=%f\n',param.alignStretchBP,param.alignStretchBP,param.alignStretchBP);
fprintf(fid,'@CLEAR\n');
fprintf(fid,'0   TBFACTOR\n');
fprintf(fid,'100  ALIGNED\n');
fprintf(fid,'@\n');
fprintf(fid,'ALIGN-ROTATION  3  KRC1=%f  KRC2=%f  KRC3=%f,\n',param.alignBendBP,param.alignBendBP,param.alignTwistBP);
fprintf(fid,'  KRC1TYPE=%s  KRC2TYPE=%s  KRC3TYPE=%s\n',param.krctypeBP,param.krctypeBP,param.krctypeBP);
fprintf(fid,'@CLEAR\n');
fprintf(fid,'0   TBFACTOR\n');
fprintf(fid,'100  ALIGNED\n');
fprintf(fid,'@\n');
fprintf(fid,'ALIGN-TRANSLATION  4  KTC1=%f  KTC2=%f  KTC3=%f\n',param.alignStretchHJ,param.alignStretchHJ,param.alignStretchHJ);
fprintf(fid,'@CLEAR\n');
fprintf(fid,'0   TBFACTOR\n');
fprintf(fid,'100  ALIGNED\n');
fprintf(fid,'@\n');
fprintf(fid,'ALIGN-ROTATION  4  KRC1=%f  KRC2=%f  KRC3=%f,\n',param.alignBendHJ_bb,0,param.alignTwistHJ_bb);
fprintf(fid,'  KRC1TYPE=%s  KRC2TYPE=%s  KRC3TYPE=%s\n',param.krctypeHJ,param.krctypeHJ,param.krctypeHJ);
fprintf(fid,'@CLEAR\n');
fprintf(fid,'0   TBFACTOR\n');
fprintf(fid,'100 ALIGNED\n');
fprintf(fid,'@\n');

% triad sets
fprintf(fid,'TRIADSETS\n');
fprintf(fid,'entries triadset aoption a1x a1y a1z a2x a2y a2z,\n');
fprintf(fid,'                 boption b0a1 b0a2 b0a3 b1a1 b1a2 b1a3 b2a1 b2a2 b2a3,\n');
fprintf(fid,'                 coption c1b1 c1b2 c1b3 c2b1 c2b2 c2b3\n');
for i=1:nElemHJE
    dirBH1 = elemHJE(i).dirBH1;
    dirBH2 = elemHJE(i).dirBH2;
    dirTH1 = elemHJE(i).dirTH1;
    dirTH2 = elemHJE(i).dirTH2;
    
    % C1_BH
    eid = 16*(i-1)+1;
    fprintf(fid,'%d vectors %f %f %f %f %f %f,\n',eid,dirBH1(1),dirBH1(2),dirBH1(3),dirBH2(1),dirBH2(2),dirBH2(3));
    fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',0.5*param.distBP,0.0,param.rHelix, 1.0,0.0,0.0, 0.0,1.0,0.0);
    fprintf(fid,'   vectors 1.0 0.0 0.0 0.0 1.0 0.0\n');
    
    % C2_BH
    eid = 16*(i-1)+2;
    fprintf(fid,'%d vectors %f %f %f %f %f %f,\n',eid,dirBH1(1),dirBH1(2),dirBH1(3),dirBH2(1),dirBH2(2),dirBH2(3));
    fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',-0.5*param.distBP,0.0,param.rHelix, 1.0,0.0,0.0, 0.0,1.0,0.0);
    fprintf(fid,'   vectors 1.0 0.0 0.0 0.0 1.0 0.0\n');
    
    % C1_TH
    eid = 16*(i-1)+3;
    fprintf(fid,'%d vectors %f %f %f %f %f %f,\n',eid,dirBH1(1),dirBH1(2),dirBH1(3),dirBH2(1),dirBH2(2),dirBH2(3));
    fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',-0.5*param.distBP*cos(param.meanAngleHJ),-0.5*param.distBP*sin(param.meanAngleHJ),-param.rHelix, 1.0,0.0,0.0, 0.0,1.0,0.0);
    fprintf(fid,'   vectors 1.0 0.0 0.0 0.0 1.0 0.0\n');
    
    % C2_TH
    eid = 16*(i-1)+4;
    fprintf(fid,'%d vectors %f %f %f %f %f %f,\n',eid,dirBH1(1),dirBH1(2),dirBH1(3),dirBH2(1),dirBH2(2),dirBH2(3));
    fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',0.5*param.distBP*cos(param.meanAngleHJ),0.5*param.distBP*sin(param.meanAngleHJ),-param.rHelix, 1.0,0.0,0.0, 0.0,1.0,0.0);
    fprintf(fid,'   vectors 1.0 0.0 0.0 0.0 1.0 0.0\n');
    
    
    % T1_BH
    eid = 16*(i-1)+5;
    angle = deg2rad(-param.angBP*(0.5+elemHJE(i).nBP_BLH-1));
    e1c = dirBH1;
    e2c = dirBH2;
    e3c = cross(e1c,e2c); e3c = e3c/norm(e3c);
    e1t = e1c;
    e2t = cos(angle)*e2c + sin(angle)*e3c;
    if ismember([elemHJE(i).idHJE,1],connHJEscaf,'rows')    % [BUG] Lots of true if no 'rows'
        fprintf(fid,'%d vectors %f %f %f %f %f %f,\n',eid,-e1t(1),-e1t(2),-e1t(3),e2t(1),e2t(2),e2t(3));
    else
        fprintf(fid,'%d vectors %f %f %f %f %f %f,\n',eid,e1t(1),e1t(2),e1t(3),e2t(1),e2t(2),e2t(3));
    end
    fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',0.0,0.0,0.0, 1.0,0.0,0.0, 0.0,1.0,0.0);
    fprintf(fid,'   vectors 1.0 0.0 0.0 0.0 1.0 0.0\n');
    
    % 3'-end of a BBE
    eid = 16*(i-1)+6;
    fprintf(fid,'%d vectors %f %f %f %f %f %f,\n',eid,e1t(1),e1t(2),e1t(3),e2t(1),e2t(2),e2t(3));
    fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',-0.5*param.distBP,param.rHelix*sin(0),-param.rHelix*cos(0), 1.0,0.0,0.0, 0.0,1.0,0.0);
%     fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',-0.5*param.distBP-10,param.rHelix*sin(0),-param.rHelix*cos(0), 1.0,0.0,0.0, 0.0,1.0,0.0);
%     fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',-0.5*param.distBP,param.rHelix*sin(0),-param.rHelix*cos(0), 1.0,0.0,0.0, 0.0,0.0,1.0);
    fprintf(fid,'   vectors 1.0 0.0 0.0 0.0 1.0 0.0\n');
    
    % 5'-end of a BBE
    eid = 16*(i-1)+7;
    fprintf(fid,'%d vectors %f %f %f %f %f %f,\n',eid,e1t(1),e1t(2),e1t(3),e2t(1),e2t(2),e2t(3));
    fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',-0.5*param.distBP,param.rHelix*sin(-a_groove),param.rHelix*cos(-a_groove), 1.0,0.0,0.0, 0.0,1.0,0.0);
%     fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',-0.5*param.distBP,param.rHelix*sin(-a_groove),param.rHelix*cos(-a_groove), -1/sqrt(2),0.0,-1/sqrt(2), 0.0,1.0,0.0);
    fprintf(fid,'   vectors 1.0 0.0 0.0 0.0 1.0 0.0\n');
    
    
    % T2_BH
    eid = 16*(i-1)+8;
    angle = deg2rad(param.angBP*(0.5+elemHJE(i).nBP_BRH-1));
    e1c = dirBH1;
    e2c = dirBH2;
    e3c = cross(e1c,e2c); e3c = e3c/norm(e3c);
    e1t = e1c;
    e2t = cos(angle)*e2c + sin(angle)*e3c;
    if ismember([elemHJE(i).idHJE,2],connHJEscaf,'rows')    % [BUG] Lots of true if no 'rows'
        fprintf(fid,'%d vectors %f %f %f %f %f %f,\n',eid,-e1t(1),-e1t(2),-e1t(3),e2t(1),e2t(2),e2t(3));
    else
        fprintf(fid,'%d vectors %f %f %f %f %f %f,\n',eid,e1t(1),e1t(2),e1t(3),e2t(1),e2t(2),e2t(3));
    end
    fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',0.0,0.0,0.0, 1.0,0.0,0.0, 0.0,1.0,0.0);
    fprintf(fid,'   vectors 1.0 0.0 0.0 0.0 1.0 0.0\n');
    
    % 3'-end of a BBE
    eid = 16*(i-1)+9;
    fprintf(fid,'%d vectors %f %f %f %f %f %f,\n',eid,e1t(1),e1t(2),e1t(3),e2t(1),e2t(2),e2t(3));
    fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',0.5*param.distBP,param.rHelix*sin(0),param.rHelix*cos(0), 1.0,0.0,0.0, 0.0,1.0,0.0);
%     fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',0.5*param.distBP+10,param.rHelix*sin(0),param.rHelix*cos(0), 1.0,0.0,0.0, 0.0,1.0,0.0);
%     fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',0.5*param.distBP,param.rHelix*sin(0),param.rHelix*cos(0), 1.0,0.0,0.0, 0.0,0.0,-1.0);
    fprintf(fid,'   vectors 1.0 0.0 0.0 0.0 1.0 0.0\n');
    
    % 5'-end of a BBE
    eid = 16*(i-1)+10;
    fprintf(fid,'%d vectors %f %f %f %f %f %f,\n',eid,e1t(1),e1t(2),e1t(3),e2t(1),e2t(2),e2t(3));
    fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',0.5*param.distBP,param.rHelix*sin(a_groove),-param.rHelix*cos(a_groove), 1.0,0.0,0.0, 0.0,1.0,0.0);
%     fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',0.5*param.distBP,param.rHelix*sin(a_groove),-param.rHelix*cos(a_groove), -1/sqrt(2),0.0,-1/sqrt(2), 0.0,1.0,0.0);
    fprintf(fid,'   vectors 1.0 0.0 0.0 0.0 1.0 0.0\n');
    
    
    % T1_TH
    eid = 16*(i-1)+11;
    angle = deg2rad(-param.angBP*(0.5+elemHJE(i).nBP_TRH-1));
    e1c = dirTH1;
    e2c = -dirTH2;
    e3c = cross(e1c,e2c); e3c = e3c/norm(e3c);
    e1t = e1c;
    e2t = cos(angle)*e2c + sin(angle)*e3c;
    if ismember([elemHJE(i).idHJE,3],connHJEscaf,'rows')    % [BUG] Lots of true if no 'rows'
        fprintf(fid,'%d vectors %f %f %f %f %f %f,\n',eid,-e1t(1),-e1t(2),-e1t(3),e2t(1),e2t(2),e2t(3));
    else
        fprintf(fid,'%d vectors %f %f %f %f %f %f,\n',eid,e1t(1),e1t(2),e1t(3),e2t(1),e2t(2),e2t(3));
    end
    fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',0.0,0.0,0.0, 1.0,0.0,0.0, 0.0,1.0,0.0);
    fprintf(fid,'   vectors 1.0 0.0 0.0 0.0 1.0 0.0\n');
    
    % 3'-end of a BBE
    eid = 16*(i-1)+12;
    fprintf(fid,'%d vectors %f %f %f %f %f %f,\n',eid,e1t(1),e1t(2),e1t(3),e2t(1),e2t(2),e2t(3));
    fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',-0.5*param.distBP,param.rHelix*sin(0),-param.rHelix*cos(0), 1.0,0.0,0.0, 0.0,1.0,0.0);
%     fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',-0.5*param.distBP-10,param.rHelix*sin(0),-param.rHelix*cos(0), 1.0,0.0,0.0, 0.0,1.0,0.0);
%     fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',-0.5*param.distBP,param.rHelix*sin(0),-param.rHelix*cos(0), 1.0,0.0,0.0, 0.0,0.0,1.0);
    fprintf(fid,'   vectors 1.0 0.0 0.0 0.0 1.0 0.0\n');
    
    % 5'-end of a BBE
    eid = 16*(i-1)+13;
    fprintf(fid,'%d vectors %f %f %f %f %f %f,\n',eid,e1t(1),e1t(2),e1t(3),e2t(1),e2t(2),e2t(3));
    fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',-0.5*param.distBP,param.rHelix*sin(-a_groove),param.rHelix*cos(-a_groove), 1.0,0.0,0.0, 0.0,1.0,0.0);
%     fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',-0.5*param.distBP,param.rHelix*sin(-a_groove),param.rHelix*cos(-a_groove), -1/sqrt(2),0.0,-1/sqrt(2), 0.0,1.0,0.0);
    fprintf(fid,'   vectors 1.0 0.0 0.0 0.0 1.0 0.0\n');
    
    
    % T2_TH
    eid = 16*(i-1)+14;
    angle = deg2rad(param.angBP*(0.5+elemHJE(i).nBP_TLH-1));
    e1c = dirTH1;
    e2c = -dirTH2;
    e3c = cross(e1c,e2c); e3c = e3c/norm(e3c);
    e1t = e1c;
    e2t = cos(angle)*e2c + sin(angle)*e3c;
    if ismember([elemHJE(i).idHJE,4],connHJEscaf,'rows')    % [BUG] Lots of true if no 'rows'
        fprintf(fid,'%d vectors %f %f %f %f %f %f,\n',eid,-e1t(1),-e1t(2),-e1t(3),e2t(1),e2t(2),e2t(3));
    else
        fprintf(fid,'%d vectors %f %f %f %f %f %f,\n',eid,e1t(1),e1t(2),e1t(3),e2t(1),e2t(2),e2t(3));
    end
    fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',0.0,0.0,0.0, 1.0,0.0,0.0, 0.0,1.0,0.0);
    fprintf(fid,'   vectors 1.0 0.0 0.0 0.0 1.0 0.0\n');
    
    % 3'-end of a BBE
    eid = 16*(i-1)+15;
    fprintf(fid,'%d vectors %f %f %f %f %f %f,\n',eid,e1t(1),e1t(2),e1t(3),e2t(1),e2t(2),e2t(3));
    fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',0.5*param.distBP,param.rHelix*sin(0),param.rHelix*cos(0), 1.0,0.0,0.0, 0.0,1.0,0.0);
%     fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',0.5*param.distBP+10,param.rHelix*sin(0),param.rHelix*cos(0), 1.0,0.0,0.0, 0.0,1.0,0.0);
%     fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',0.5*param.distBP,param.rHelix*sin(0),param.rHelix*cos(0), 1.0,0.0,0.0, 0.0,0.0,-1.0);
    fprintf(fid,'   vectors 1.0 0.0 0.0 0.0 1.0 0.0\n');
    
    % 5'-end of a BBE
    eid = 16*(i-1)+16;
    fprintf(fid,'%d vectors %f %f %f %f %f %f,\n',eid,e1t(1),e1t(2),e1t(3),e2t(1),e2t(2),e2t(3));
    fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',0.5*param.distBP,param.rHelix*sin(a_groove),-param.rHelix*cos(a_groove), 1.0,0.0,0.0, 0.0,1.0,0.0);
%     fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',0.5*param.distBP,param.rHelix*sin(a_groove),-param.rHelix*cos(a_groove), -1/sqrt(2),0.0,-1/sqrt(2), 0.0,1.0,0.0);
    fprintf(fid,'   vectors 1.0 0.0 0.0 0.0 1.0 0.0\n');
end
for i=1:nElemDSE
    dirBH1 = elemDSE(i).dirBH1;
    dirBH2 = elemDSE(i).dirBH2;
    
    % T1_BH
    eid = 16*nElemHJE + 6*(i-1) + 1;
    angle = 0;
    e1c = dirBH1;
    e2c = dirBH2;
    e3c = cross(e1c,e2c); e3c = e3c/norm(e3c);
    e1t = e1c;
    e2t = cos(angle)*e2c + sin(angle)*e3c;
    if ismember([elemDSE(i).idHJE,1],connHJEscaf,'rows')
        fprintf(fid,'%d vectors %f %f %f %f %f %f,\n',eid,-e1t(1),-e1t(2),-e1t(3),e2t(1),e2t(2),e2t(3));
    else
        fprintf(fid,'%d vectors %f %f %f %f %f %f,\n',eid,e1t(1),e1t(2),e1t(3),e2t(1),e2t(2),e2t(3));
    end
    fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',0.0,0.0,0.0, 1.0,0.0,0.0, 0.0,1.0,0.0);
    fprintf(fid,'   vectors 1.0 0.0 0.0 0.0 1.0 0.0\n');
    
    % 3'-end of a BBE
    eid = 16*nElemHJE + 6*(i-1) + 2;
    fprintf(fid,'%d vectors %f %f %f %f %f %f,\n',eid,e1t(1),e1t(2),e1t(3),e2t(1),e2t(2),e2t(3));
    fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',-0.5*param.distBP,param.rHelix*sin(0),-param.rHelix*cos(0), 1.0,0.0,0.0, 0.0,1.0,0.0);
%     fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',-0.5*param.distBP,0.0,-param.rHelix, 1.0,0.0,0.0, 0.0,1.0,0.0);
    fprintf(fid,'   vectors 1.0 0.0 0.0 0.0 1.0 0.0\n');
    
    % 5'-end of a BBE
    eid = 16*nElemHJE + 6*(i-1) + 3;
    fprintf(fid,'%d vectors %f %f %f %f %f %f,\n',eid,e1t(1),e1t(2),e1t(3),e2t(1),e2t(2),e2t(3));
    fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',-0.5*param.distBP,param.rHelix*sin(-a_groove),param.rHelix*cos(-a_groove), 1.0,0.0,0.0, 0.0,1.0,0.0);
%     fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',-0.5*param.distBP,0.0,param.rHelix, 1.0,0.0,0.0, 0.0,1.0,0.0);
    fprintf(fid,'   vectors 1.0 0.0 0.0 0.0 1.0 0.0\n');
    
    
    % T2_BH
    eid = 16*nElemHJE + 6*(i-1) + 4;
    angle = deg2rad(param.angBP*(elemDSE(i).nBP-1));
    e1c = dirBH1;
    e2c = dirBH2;
    e3c = cross(e1c,e2c); e3c = e3c/norm(e3c);
    e1t = e1c;
    e2t = cos(angle)*e2c + sin(angle)*e3c;
    if ismember([elemDSE(i).idHJE,2],connHJEscaf,'rows')
        fprintf(fid,'%d vectors %f %f %f %f %f %f,\n',eid,-e1t(1),-e1t(2),-e1t(3),e2t(1),e2t(2),e2t(3));
    else
        fprintf(fid,'%d vectors %f %f %f %f %f %f,\n',eid,e1t(1),e1t(2),e1t(3),e2t(1),e2t(2),e2t(3));
    end
    fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',0.0,0.0,0.0, 1.0,0.0,0.0, 0.0,1.0,0.0);
    fprintf(fid,'   vectors 1.0 0.0 0.0 0.0 1.0 0.0\n');
    
    % 3'-end of a BBE
    eid = 16*nElemHJE + 6*(i-1) + 5;
    fprintf(fid,'%d vectors %f %f %f %f %f %f,\n',eid,e1t(1),e1t(2),e1t(3),e2t(1),e2t(2),e2t(3));
    fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',0.5*param.distBP,param.rHelix*sin(0),param.rHelix*cos(0), 1.0,0.0,0.0, 0.0,1.0,0.0);
%     fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',0.5*param.distBP,0.0,param.rHelix, 1.0,0.0,0.0, 0.0,1.0,0.0);
    fprintf(fid,'   vectors 1.0 0.0 0.0 0.0 1.0 0.0\n');
    
    % 5'-end of a BBE
    eid = 16*nElemHJE + 6*(i-1) + 6;
    fprintf(fid,'%d vectors %f %f %f %f %f %f,\n',eid,e1t(1),e1t(2),e1t(3),e2t(1),e2t(2),e2t(3));
    fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',0.5*param.distBP,param.rHelix*sin(a_groove),-param.rHelix*cos(a_groove), 1.0,0.0,0.0, 0.0,1.0,0.0);
%     fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',0.5*param.distBP,0.0,-param.rHelix, 1.0,0.0,0.0, 0.0,1.0,0.0);
    fprintf(fid,'   vectors 1.0 0.0 0.0 0.0 1.0 0.0\n');
end

% open nicks
FEModel.nickOpen = [];
eid = 16*nElemHJE + 6*nElemDSE;

for i=1:nElemHJE
    dirBH1 = elemHJE(i).dirBH1;
    dirBH2 = elemHJE(i).dirBH2;
    dirTH1 = elemHJE(i).dirTH1;
    dirTH2 = elemHJE(i).dirTH2;
    
    angle = deg2rad(-param.angBP*(0.5+elemHJE(i).nBP_BLH-1));
    e1c = dirBH1;
    e2c = dirBH2;
    e3c = cross(e1c,e2c); e3c = e3c/norm(e3c);
    e1t_BH = e1c;
    e2t_BH = cos(angle)*e2c + sin(angle)*e3c;
    e3t_BH = cross(e1t_BH,e2t_BH); e3t_BH = e3t_BH/norm(e3t_BH);
    
    angle = deg2rad(-param.angBP*(0.5+elemHJE(i).nBP_TRH-1));
    e1c = dirTH1;
    e2c = -dirTH2;
    e3c = cross(e1c,e2c); e3c = e3c/norm(e3c);
    e1t_TH = e1c;
    e2t_TH = cos(angle)*e2c + sin(angle)*e3c;
    e3t_TH = cross(e1t_TH,e2t_TH); e3t_TH = e3t_TH/norm(e3t_TH);
    
    nickOpen = nick(nick(:,1)==elemHJE(i).idHJE & nick(:,3)>0, 2:3);
    
    for j = 1 : size(nickOpen,1)
        nid = elemHJE(i).nidT1_BH - 1 + nickOpen(j,1);
        
        if(nickOpen(j,1) <= elemHJE(i).nBP_BLH + elemHJE(i).nBP_BRH);
            j1 = nickOpen(j,1);
            angle = deg2rad(param.angBP * (0.5 + j1 - 1));
            e1n = e1t_BH;
            e2n = cos(angle)*e2t_BH + sin(angle)*e3t_BH;
        else
            j1 = nickOpen(j,1) - (elemHJE(i).nBP_BLH + elemHJE(i).nBP_BRH);
            angle = deg2rad(param.angBP * (0.5 + j1 - 1));
            e1n = e1t_TH;
            e2n = cos(angle)*e2t_TH + sin(angle)*e3t_TH;
        end
        
        eid = eid + 1;
        fprintf(fid,'%d vectors %f %f %f %f %f %f,\n',eid,e1n(1),e1n(2),e1n(3),e2n(1),e2n(2),e2n(3));
        if(nickOpen(j,2) == 1)
            fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',0.5*param.distBP,0.0,-param.rHelix, 1.0,0.0,0.0, 0.0,1.0,0.0);
        elseif(nickOpen(j,2) == 2)
            fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',0.5*param.distBP,0.0,param.rHelix, 1.0,0.0,0.0, 0.0,1.0,0.0);
        else
            error('Exception.');
        end
        fprintf(fid,'   vectors 1.0 0.0 0.0 0.0 1.0 0.0\n');
        
        eid = eid + 1;
        fprintf(fid,'%d vectors %f %f %f %f %f %f,\n',eid,e1n(1),e1n(2),e1n(3),e2n(1),e2n(2),e2n(3));
        if(nickOpen(j,2) == 1)
            fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',-0.5*param.distBP,0.0,-param.rHelix, 1.0,0.0,0.0, 0.0,1.0,0.0);
        elseif(nickOpen(j,2) == 2)
            fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',-0.5*param.distBP,0.0,param.rHelix, 1.0,0.0,0.0, 0.0,1.0,0.0);
        else
            error('Exception.');
        end
        fprintf(fid,'   vectors 1.0 0.0 0.0 0.0 1.0 0.0\n');
        
        FEModel.nickOpen = cat(1, FEModel.nickOpen, [nid, eid-1, nid+1, eid]);
    end
end

for i=1:nElemDSE
    dirBH1 = elemDSE(i).dirBH1;
    dirBH2 = elemDSE(i).dirBH2;
    
    angle = 0;
    e1c = dirBH1;
    e2c = dirBH2;
    e3c = cross(e1c,e2c); e3c = e3c/norm(e3c);
    e1t = e1c;
    e2t = cos(angle)*e2c + sin(angle)*e3c;
    e3t = cross(e1t,e2t); e3t = e3t/norm(e3t);
    
    nickOpen = nick(nick(:,1)==elemDSE(i).idHJE & nick(:,3)>0, 2:3);
    
    for j = 1 : size(nickOpen,1)
        nid = elemDSE(i).nidT1 - 1 + nickOpen(j,1);
        
        j1 = nickOpen(j,1);
        angle = deg2rad(param.angBP * (0.5 + j1 - 1));
        e1n = e1t;
        e2n = cos(angle)*e2t + sin(angle)*e3t;
        
        eid = eid + 1;
        fprintf(fid,'%d vectors %f %f %f %f %f %f,\n',eid,e1n(1),e1n(2),e1n(3),e2n(1),e2n(2),e2n(3));
        if(nickOpen(j,2) == 1)
            fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',0.5*param.distBP,0.0,-param.rHelix, 1.0,0.0,0.0, 0.0,1.0,0.0);
        elseif(nickOpen(j,2) == 2)
            fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',0.5*param.distBP,0.0,param.rHelix, 1.0,0.0,0.0, 0.0,1.0,0.0);
        else
            error('Exception.');
        end
        fprintf(fid,'   vectors 1.0 0.0 0.0 0.0 1.0 0.0\n');
        
        eid = eid + 1;
        fprintf(fid,'%d vectors %f %f %f %f %f %f,\n',eid,e1n(1),e1n(2),e1n(3),e2n(1),e2n(2),e2n(3));
        if(nickOpen(j,2) == 1)
            fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',-0.5*param.distBP,0.0,-param.rHelix, 1.0,0.0,0.0, 0.0,1.0,0.0);
        elseif(nickOpen(j,2) == 2)
            fprintf(fid,'   vectors %f %f %f %f %f %f %f %f %f,\n',-0.5*param.distBP,0.0,param.rHelix, 1.0,0.0,0.0, 0.0,1.0,0.0);
        else
            error('Exception.');
        end
        fprintf(fid,'   vectors 1.0 0.0 0.0 0.0 1.0 0.0\n');
        
        FEModel.nickOpen = cat(1, FEModel.nickOpen, [nid, eid-1, nid+1, eid]);
    end
end


% set triad sets to nodes
fprintf(fid,'SET-TRIADSET-NODES\n');
fprintf(fid,'@CLEAR\n');
for i=1:nElemHJE
    
    % C1_BH
    eid = 16*(i-1)+1;
    nid = elemHJE(i).nidC1_BH;
    fprintf(fid,'%d  %d\n',nid,eid);
    
    % C2_BH
    eid = 16*(i-1)+2;
    nid = elemHJE(i).nidC2_BH;
    fprintf(fid,'%d  %d\n',nid,eid);
    
    % C1_TH
    eid = 16*(i-1)+3;
    nid = elemHJE(i).nidC1_TH;
    fprintf(fid,'%d  %d\n',nid,eid);
    
    % C2_TH
    eid = 16*(i-1)+4;
    nid = elemHJE(i).nidC2_TH;
    fprintf(fid,'%d  %d\n',nid,eid);
    
    % T1_BH
    eid = 16*(i-1)+5;
    nid = elemHJE(i).nidT1_BH;
    fprintf(fid,'%d  %d\n',nid,eid);
    
    % T2_BH
    eid = 16*(i-1)+8;
    nid = elemHJE(i).nidT2_BH;
    fprintf(fid,'%d  %d\n',nid,eid);
    
    % T1_TH
    eid = 16*(i-1)+11;
    nid = elemHJE(i).nidT1_TH;
    fprintf(fid,'%d  %d\n',nid,eid);
    
    % T2_TH
    eid = 16*(i-1)+14;
    nid = elemHJE(i).nidT2_TH;
    fprintf(fid,'%d  %d\n',nid,eid);
    
end
for i=1:nElemDSE
    
    % T1_BH
    eid = 16*nElemHJE + 6*(i-1)+1;
    nid = elemDSE(i).nidT1;
    fprintf(fid,'%d  %d\n',nid,eid);
    
    % T2_BH
    eid = 16*nElemHJE + 6*(i-1)+4;
    nid = elemDSE(i).nidT2;
    fprintf(fid,'%d  %d\n',nid,eid);
    
end
for i = 1 : size(FEModel.BBE, 1)
    fprintf(fid,'%d  %d\n',FEModel.BBE(i,1),FEModel.BBE(i,2));
    fprintf(fid,'%d  %d\n',FEModel.BBE(i,3),FEModel.BBE(i,4));
end
for i = 1 : size(FEModel.nickOpen, 1)
    fprintf(fid,'%d  %d\n',FEModel.nickOpen(i,1),FEModel.nickOpen(i,2));
    fprintf(fid,'%d  %d\n',FEModel.nickOpen(i,3),FEModel.nickOpen(i,4));
end
fprintf(fid,'@\n');

% alignment elements for 4-way junctions
FEModel.alignHJE = zeros(2*nElemHJE,2);
FEModel.typeHJE = cell(2*nElemHJE,1);

% alignment elements for double-stranded crossovers
iElemHJE_ds = 0;
if(nElemHJE_ds > 0)
    fprintf(fid,'EGROUP ALIGNMENT NAME=%d ALIGN-TRANSLATION=1 ALIGN-ROTATION=1 RESULTS=ALIGNMENTS\n',4+nEntSprings+1);
    fprintf(fid,'ENODES GROUP=%d\n',4+nEntSprings+1);
    fprintf(fid,'@CLEAR\n');
%     iElemHJE_ds = 0;
    for i=1:nElemHJE
        if(strcmp(elemHJE(i).type, 'ds'))
            fprintf(fid,'%d  %d  %d\n',nElemDSDNA+nEntSprings+2*iElemHJE_ds+1,elemHJE(i).nidC1_BH,elemHJE(i).nidC2_TH);
            fprintf(fid,'%d  %d  %d\n',nElemDSDNA+nEntSprings+2*iElemHJE_ds+2,elemHJE(i).nidC2_BH,elemHJE(i).nidC1_TH);
            FEModel.alignHJE(2*(i-1)+1,:) = [elemHJE(i).nidC1_BH, elemHJE(i).nidC2_TH];
            FEModel.alignHJE(2*(i-1)+2,:) = [elemHJE(i).nidC2_BH, elemHJE(i).nidC1_TH];
            FEModel.typeHJE{2*(i-1)+1} = 'ds';
            FEModel.typeHJE{2*(i-1)+2} = 'ds';
            iElemHJE_ds = iElemHJE_ds + 1;
        end
    end
    fprintf(fid,'@\n');
    fprintf(fid,'EDATA GROUP=%d\n',4+nEntSprings+1);
    % elem #, print flag, save flag, death time, align translation, align distance, align rotation, triadset 1, triadset 2
    fprintf(fid,'@CLEAR\n');
    iElemHJE_ds = 0;
    for i=1:nElemHJE
        if(strcmp(elemHJE(i).type, 'ds'))
            fprintf(fid,'%d  ''Default''  ''Default''  0.0  1  0  1  %d  %d\n',nElemDSDNA+nEntSprings+2*iElemHJE_ds+1,16*(i-1)+1,16*(i-1)+4);
            fprintf(fid,'%d  ''Default''  ''Default''  0.0  1  0  1  %d  %d\n',nElemDSDNA+nEntSprings+2*iElemHJE_ds+2,16*(i-1)+2,16*(i-1)+3);
            iElemHJE_ds = iElemHJE_ds + 1;
        end
    end
    fprintf(fid,'@\n');
end

% alignment elements for single-stranded crossovers
if(nElemHJE_ss > 0)
    fprintf(fid,'EGROUP ALIGNMENT NAME=%d ALIGN-TRANSLATION=2 ALIGN-ROTATION=2 RESULTS=ALIGNMENTS\n',5+nEntSprings+1);
    fprintf(fid,'ENODES GROUP=%d\n',5+nEntSprings+1);
    fprintf(fid,'@CLEAR\n');
    iElemHJE_ss = iElemHJE_ds;
    for i=1:nElemHJE
        if(strcmp(elemHJE(i).type, 'ss'))
            fprintf(fid,'%d  %d  %d\n',nElemDSDNA+nEntSprings+2*iElemHJE_ss+1,elemHJE(i).nidC1_BH,elemHJE(i).nidC2_TH);
            fprintf(fid,'%d  %d  %d\n',nElemDSDNA+nEntSprings+2*iElemHJE_ss+2,elemHJE(i).nidC2_BH,elemHJE(i).nidC1_TH);
            FEModel.alignHJE(2*(i-1)+1,:) = [elemHJE(i).nidC1_BH, elemHJE(i).nidC2_TH];
            FEModel.alignHJE(2*(i-1)+2,:) = [elemHJE(i).nidC2_BH, elemHJE(i).nidC1_TH];
            FEModel.typeHJE{2*(i-1)+1} = 'ss';
            FEModel.typeHJE{2*(i-1)+2} = 'ss';
            iElemHJE_ss = iElemHJE_ss + 1;
        end
    end
    fprintf(fid,'@\n');
    fprintf(fid,'EDATA GROUP=%d\n',5+nEntSprings+1);
    % elem #, print flag, save flag, death time, align translation, align distance, align rotation, triadset 1, triadset 2
    fprintf(fid,'@CLEAR\n');
    iElemHJE_ss = iElemHJE_ds;
    for i=1:nElemHJE
        if(strcmp(elemHJE(i).type, 'ss'))
            fprintf(fid,'%d  ''Default''  ''Default''  0.0  2  0  2  %d  %d\n',nElemDSDNA+nEntSprings+2*iElemHJE_ss+1,16*(i-1)+1,16*(i-1)+4);
            fprintf(fid,'%d  ''Default''  ''Default''  0.0  2  0  2  %d  %d\n',nElemDSDNA+nEntSprings+2*iElemHJE_ss+2,16*(i-1)+2,16*(i-1)+3);
            iElemHJE_ss = iElemHJE_ss + 1;
        end
    end
    fprintf(fid,'@\n');
end

% alignment elements for HJ element unit connection
if size(connHJE,1)==0
    FEModel.alignDSDNA = [];
else
    FEModel.alignDSDNA = zeros(size(connHJE,1),2);
    fprintf(fid,'EGROUP ALIGNMENT NAME=%d ALIGN-TRANSLATION=3 ALIGN-ROTATION=3 RESULTS=ALIGNMENTS\n',6+nEntSprings+1);
    fprintf(fid,'ENODES GROUP=%d\n',6+nEntSprings+1);
    fprintf(fid,'@CLEAR\n');
    for i=1:size(connHJE,1)
        
        hjeid1 = connHJE(i,1);
        if HJE(hjeid1).nWay == 4
            idx1 = find([elemHJE.idHJE]==hjeid1,1);
            switch(connHJE(i,2))
                case 1,
                    nid1 = elemHJE(idx1).nidT1_BH;
                case 2,
                    nid1 = elemHJE(idx1).nidT2_BH;
                case 3,
                    nid1 = elemHJE(idx1).nidT1_TH;
                case 4,
                    nid1 = elemHJE(idx1).nidT2_TH;
            end
        elseif HJE(hjeid1).nWay == 2
            idx1 = find([elemDSE.idHJE]==hjeid1,1);
            switch(connHJE(i,2))
                case 1,
                    nid1 = elemDSE(idx1).nidT1;
                case 2,
                    nid1 = elemDSE(idx1).nidT2;
            end
        else
            error('Unknown %d-way junction!\n',HJE(hjeid1).nWay);
        end
        
        hjeid2 = connHJE(i,3);
        if HJE(hjeid2).nWay == 4
            idx2 = find([elemHJE.idHJE]==hjeid2,1);
            switch(connHJE(i,4))
                case 1,
                    nid2 = elemHJE(idx2).nidT1_BH;
                case 2,
                    nid2 = elemHJE(idx2).nidT2_BH;
                case 3,
                    nid2 = elemHJE(idx2).nidT1_TH;
                case 4,
                    nid2 = elemHJE(idx2).nidT2_TH;
            end
        elseif HJE(hjeid2).nWay == 2
            idx2 = find([elemDSE.idHJE]==hjeid2,1);
            switch(connHJE(i,4))
                case 1,
                    nid1 = elemDSE(idx2).nidT1;
                case 2,
                    nid1 = elemDSE(idx2).nidT2;
            end
        else
            error('Unknown %d-way junction!\n',HJE(hjeid2).nWay);
        end
        
        fprintf(fid,'%d  %d  %d\n',nElemDSDNA+nEntSprings+2*nElemHJE+i,nid1,nid2);
        FEModel.alignDSDNA(i,:) = [nid1, nid2];
    end
    fprintf(fid,'@\n');
end

if size(connHJE,1)~=0
    fprintf(fid,'EDATA GROUP=%d\n',6+nEntSprings+1);
    fprintf(fid,'@CLEAR\n');
    for i=1:size(connHJE,1)
        
        hjeid1 = connHJE(i,1);
        if HJE(hjeid1).nWay == 4
            idx1 = find([elemHJE.idHJE]==hjeid1,1);
            switch(connHJE(i,2))
                case 1,
                    eid1 = 16*(idx1-1)+5;
                case 2,
                    eid1 = 16*(idx1-1)+8;
                case 3,
                    eid1 = 16*(idx1-1)+11;
                case 4,
                    eid1 = 16*(idx1-1)+14;
            end
        elseif HJE(hjeid1).nWay == 2
            idx1 = find([elemDSE.idHJE]==hjeid1,1);
            switch(connHJE(i,2))
                case 1,
                    eid1 = 16*nElemHJE + 6*(idx1-1)+1;
                case 2,
                    eid1 = 16*nElemHJE + 6*(idx1-1)+4;
            end
        else
            error('Unknown %d-way junction!\n',HJE(hjeid1).nWay);
        end
        
        hjeid2 = connHJE(i,3);
        if HJE(hjeid2).nWay == 4
            idx2 = find([elemHJE.idHJE]==hjeid2,1);
            switch(connHJE(i,4))
                case 1,
                    eid2 = 16*(idx2-1)+5;
                case 2,
                    eid2 = 16*(idx2-1)+8;
                case 3,
                    eid2 = 16*(idx2-1)+11;
                case 4,
                    eid2 = 16*(idx2-1)+14;
            end
        elseif HJE(hjeid2).nWay == 2
            idx2 = find([elemDSE.idHJE]==hjeid2,1);
            switch(connHJE(i,4))
                case 1,
                    eid2 = 16*nElemHJE + 6*(idx2-1)+1;
                case 2,
                    eid2 = 16*nElemHJE + 6*(idx2-1)+4;
            end
        else
            error('Unknown %d-way junction!\n',HJE(hjeid2).nWay);
        end
        
        fprintf(fid,'%d  ''Default''  ''Default''  0.0  3  0  3  %d  %d\n',nElemDSDNA+nEntSprings+2*nElemHJE+i,eid1,eid2);
        
    end
    fprintf(fid,'@\n');
end

% alignment elements for bulges and open nicks phosphodiester backbones
if(nBBE > 0 || size(FEModel.nickOpen,1) > 0)
    fprintf(fid,'EGROUP ALIGNMENT NAME=%d ALIGN-TRANSLATION=4 ALIGN-ROTATION=4 RESULTS=ALIGNMENTS\n',7+nEntSprings+1);
    fprintf(fid,'ENODES GROUP=%d\n',7+nEntSprings+1);
    fprintf(fid,'@CLEAR\n');
    for i = 1:nBBE
        fprintf(fid,'%d  %d  %d\n', nElemDSDNA + nEntSprings + 2*nElemHJE_ds + 2*nElemHJE_ss + size(connHJE,1) + i, FEModel.BBE(i,1), FEModel.BBE(i,3));
    end
    for i = 1:size(FEModel.nickOpen,1)
        fprintf(fid,'%d  %d  %d\n', nElemDSDNA + nEntSprings + 2*nElemHJE_ds + 2*nElemHJE_ss + size(connHJE,1) + nBBE + i, FEModel.nickOpen(i,1), FEModel.nickOpen(i,3));
    end
    fprintf(fid,'@\n');
    fprintf(fid,'EDATA GROUP=%d\n',7+nEntSprings+1);
    % elem #, print flag, save flag, death time, align translation, align distance, align rotation, triadset 1, triadset 2
    fprintf(fid,'@CLEAR\n');
    for i = 1:nBBE
        fprintf(fid,'%d  ''Default''  ''Default''  0.0  4  0  4  %d  %d\n', nElemDSDNA + nEntSprings + 2*nElemHJE_ds + 2*nElemHJE_ss + size(connHJE,1) + i, FEModel.BBE(i,2), FEModel.BBE(i,4));
    end
    for i = 1:size(FEModel.nickOpen,1)
        fprintf(fid,'%d  ''Default''  ''Default''  0.0  4  0  4  %d  %d\n', nElemDSDNA + nEntSprings + 2*nElemHJE_ds + 2*nElemHJE_ss + size(connHJE,1) + nBBE + i, FEModel.nickOpen(i,2), FEModel.nickOpen(i,4));
    end
    fprintf(fid,'@\n');
end

fprintf(fid,'MASS-MATRIX LUMPED ETA=1.0\n');

fprintf(fid,'ITERATION METHOD=FULL-NEWTON LINE-SEA=DEFAULT MAX-ITER=20,\n');
fprintf(fid,'     PRINTOUT=NONE  PLASTIC-=1\n');
fprintf(fid,'TOLERANCES ITERATION CONVERGE=EF ETOL=1.00000000000000E-06,\n');
fprintf(fid,'     RTOL=0.0100000000000000 RNORM=10.000000000000,\n');
fprintf(fid,'     RMNORM=100.00000000000 RCTOL=0.0500000000000000,\n');
fprintf(fid,'     STOL=0.500000000000000 RCONSM=0.0100000000000000,\n');
fprintf(fid,'     ENLSTH=0.00000000000000 LSLOWER=0.00100000000000000,\n');
fprintf(fid,'     LSUPPER=0.00000000000000 MAXDISP=0.00000000000000\n');
fprintf(fid,'PRINTOUT ECHO=NO PRINTDEF=YES INPUT-DA=4 OUTPUT=SELECTED DISPLACE=YES,\n');
fprintf(fid,'     VELOCITI=YES ACCELERA=YES IDISP=NO ITEMP=NO ISTRAIN=NO IPIPE=NO,\n');
fprintf(fid,'     STORAGE=NO LARGE-ST=NONE ENERGIES=YES\n');
fprintf(fid,'NODESAVE-STE ELEMSAVE=OVERWRITE\n');
fprintf(fid,'@CLEAR\n');
for i=1:size(param.adinaSaveTimeStep,1)
    fprintf(fid,'%d  %d  %d  %d\n',param.adinaSaveTimeStep(i,1),param.adinaSaveTimeStep(i,2),param.adinaSaveTimeStep(i,3),param.adinaSaveTimeStep(i,4));
end
fprintf(fid,'@\n');
fprintf(fid,'PRINT-STEPS SUBSTRUCT=0 REUSE=1\n');
fprintf(fid,'@CLEAR\n');
for i=1:size(param.adinaSaveTimeStep,1)
    fprintf(fid,'%d  %d  %d  %d\n',param.adinaSaveTimeStep(i,1),param.adinaSaveTimeStep(i,2),param.adinaSaveTimeStep(i,3),param.adinaSaveTimeStep(i,4));
end
fprintf(fid,'@\n');

fprintf(fid,'MASTER ANALYSIS=DYNAMIC MODEX=EXECUTE TSTART=0.00000000000000 IDOF=0,\n');
fprintf(fid,'     OVALIZAT=NONE FLUIDPOT=AUTOMATIC CYCLICPA=1 IPOSIT=CONTINUE,\n');
fprintf(fid,'     REACTION=YES INITIALS=NO FSINTERA=NO IRINT=DEFAULT CMASS=NO,\n');
fprintf(fid,'     SHELLNDO=AUTOMATIC AUTOMATI=ATS SOLVER=NONSYM-SP,\n');
fprintf(fid,'     CONTACT-=CONSTRAINT-FUNCTION TRELEASE=0.00000000000000,\n');
fprintf(fid,'     RESTART-=NO FRACTURE=NO LOAD-CAS=NO LOAD-PEN=NO SINGULAR=YES,\n');
fprintf(fid,'     STIFFNES=0.000100000000000000 MAP-OUTP=NONE MAP-FORM=NO,\n');
fprintf(fid,'     NODAL-DE='''' POROUS-C=NO ADAPTIVE=0 ZOOM-LAB=1 AXIS-CYC=0,\n');
fprintf(fid,'     PERIODIC=NO VECTOR-S=GEOMETRY EPSI-FIR=NO STABILIZ=NO,\n');
fprintf(fid,'     STABFACT=1.00000000000000E-10 RESULTS=PORTHOLE FEFCORR=NO,\n');
fprintf(fid,'     BOLTSTEP=1 EXTEND-S=YES CONVERT-=NO DEGEN=YES TMC-MODE=NO,\n');
fprintf(fid,'     ENSIGHT-=NO\n');
fprintf(fid,'AUTOMATIC TIME-STEPPING MAXSUBD=1000000000\n');
fprintf(fid,'ANALYSIS DYNAMIC METHOD=BATHE\n');
fprintf(fid,'RAYLEIGH-DAMPING alpha=%f beta=%f UPDATE=NONE\n',param.dampingMass,param.dampingStiff);
fprintf(fid,'DATAEND\n');

% save files
fprintf(fid,'DATABASE SAVE PERMFILE=,\n');
fprintf(fid,'''');
fprintf(fid,'%s',fullfile(tarDIR,strcat(bodyFN,'.idb')));
fprintf(fid,''',\n');
fprintf(fid,'     PROMPT=NO\n');
fprintf(fid,'ADINA OPTIMIZE=SOLVER FILE=,\n');
fprintf(fid,'''');
fprintf(fid,'%s',fullfile(tarDIR,strcat(bodyFN,'.dat')));
fprintf(fid,''',\n');
fprintf(fid,'     FIXBOUND=YES MIDNODE=NO OVERWRIT=YES\n');
fprintf(fid,'EXIT SAVE=NO IMMEDIATE=YES\n');

fclose(fid);

FEModel.nodalTriad = nodalTriad;
save(fullfile(tarDIR,strcat(bodyFN,'_FEModel.mat')),'FEModel','param','elemHJE','elemDSE','connHJE','connSSE');
%-------------------------------------------------------------------


% ADINA: .in --> .dat
fprintf(1,'\t\t+ Translating: .in --> .dat\n');
adinaIN = fullfile(tarDIR,strcat(bodyFN,'.in'));
runAUI = sprintf('%s %s %s',param.auiEXE,param.auiOPTION,adinaIN);
system(runAUI);

% ADINA: .dat --> .por
fprintf(1,'\t\t+ Solving\n');
adinaDAT = fullfile(tarDIR,strcat(bodyFN,'.dat'));
runADINA = sprintf('%s %s %s',param.adinaEXE,param.adinaOPTION,adinaDAT);
system(runADINA);

end
