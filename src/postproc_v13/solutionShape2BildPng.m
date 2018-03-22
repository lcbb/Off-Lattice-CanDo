% Five new files will be created in the working directory:
%     - *.bild
%     - *_chimeraScr.py
%     - *_view<1,2,3>.png
%
function [] = solutionShape2BildPng(matPATH,pNode,suffixFN,nodeVal,elemVal,valRange,triad,cmap,param,outDIR)

rTruss = 2.0;   % Unit: Angstrom
lenTriad = 20;  % Unit: Angstrom
rTriad = 1;     % Unit: Angstrom

% Read in the finite element model
load(matPATH,'FEModel');

% Get the prefix of the filenames
[tarDIR,bodyFN] = fileparts(matPATH);
bodyFN = bodyFN(1:end-8);   % remove the suffix '_FEModel'
assert(~isempty(bodyFN));

% A list of beam elements to be rendered
elemRender = [ ...
               FEModel.beamsNormal; ...
               FEModel.beamsNick; ...
%                FEModel.shrinkBeamsNormal; ...
%                FEModel.shrinkBeamsNick; ...
             ];
nElemRender = size(elemRender,1);

% A list of truss elements to be rendered
if(isfield(FEModel,'rigidLinks'))
    elemRenderTruss = [ ...
                       FEModel.rigidLinks; ...
                       FEModel.entSprings; ...
                      ];
else
    elemRenderTruss = [FEModel.entSprings];
end
nElemRenderTruss = size(elemRenderTruss,1);

% Check the number of nodes and beam elements
nFrame = size(pNode,3);
assert(isempty(nodeVal) || isempty(elemVal));
if(~isempty(nodeVal))
    assert(size(nodeVal,1)==size(FEModel.node,1) && size(nodeVal,2)==nFrame);
end
if(~isempty(elemVal))
    assert(size(elemVal,1)==size(elemRender,1) && size(elemVal,2)==nFrame);
end

% Create the working directory
if(exist(fullfile(tarDIR,outDIR),'dir')~=7)
    mkdir(fullfile(tarDIR,outDIR));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computing cylinder information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cylinder = cell(nFrame,1);
for loop = 1:nFrame
    iCylinder = 0;
    for i=1:nElemRender             % ... for beam elements
        nid1 = elemRender(i,2);
        nid2 = elemRender(i,3);
        coord1 = pNode(nid1,:,loop);
        coord2 = pNode(nid2,:,loop);
        if norm(coord1-coord2)<1e-4
            continue;
        end
        coordM = 0.5*(coord1+coord2);
        iCylinder = iCylinder + 1;
        cylinder{loop}(iCylinder,:) = [nid1, coord1, coordM, param.rHelix];
        iCylinder = iCylinder + 1;
        cylinder{loop}(iCylinder,:) = [nid2, coordM, coord2, param.rHelix];
    end
    for j=1:nElemRenderTruss        % ... for truss elements
        nid1 = elemRenderTruss(j,2);
        nid2 = elemRenderTruss(j,3);
        if nid1 ~= 0 && nid2 ~= 0 && nid1~=nid2
            coord1 = pNode(nid1,:,loop);
            coord2 = pNode(nid2,:,loop);
            if norm(coord1-coord2)<1e-4
                continue;
            end
            coordM = 0.5*(coord1+coord2);
            iCylinder = iCylinder + 1;
            cylinder{loop}(iCylinder,:) = [nid1, coord1, coordM, rTruss];
            iCylinder = iCylinder + 1;
            cylinder{loop}(iCylinder,:) = [nid2, coordM, coord2, rTruss];
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computing RGB information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(~isempty(nodeVal))
    entry = nodeVal;
elseif(~isempty(elemVal))
    entry = elemVal;
else
    entry = [];
    assert(isempty(valRange));
    assert(isempty(cmap));
end

if(~isempty(valRange))
    minEntry = valRange(1);
    maxEntry = valRange(2);
else
    minEntry = min(entry(:));
    maxEntry = max(entry(:));
end


if(~isempty(entry))
    relEntry = (entry-minEntry)/(maxEntry-minEntry);
    %maxEntry = quantile(relEntry(:),0.95); % removing high peaks
    relEntry(relEntry<0) = 0;
    relEntry(relEntry>1) = 1;
    %bound = [0:0.005:maxEntry 1.0];
    bound = 0:0.005:1.0;
    if(strcmp(cmap,'jet'))
        color = jet(numel(bound));
    elseif(strcmp(cmap,'hot'))
        color = hot(numel(bound));
    elseif(strcmp(cmap,'cool'))
        color = cool(numel(bound));
    elseif(strcmp(cmap,'red'))
        color = cmapGen(0, 0.75, numel(bound), 'single');
    elseif(strcmp(cmap,'orange_blue'))
        color = cmapGen([0.575 0.075], 0.75, numel(bound), 'double');
    else
        error('Undefined color map.');
    end
    [~,BIN] = histc(relEntry,bound);
else
    BIN = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generating the BILD file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bildPATH = cell(nFrame,1);
for loop = 1:nFrame
    if(nFrame==1)
        bildPATH{loop} = fullfile(tarDIR, outDIR, strcat(bodyFN,suffixFN,'.bild'));
    else
        bildPATH{loop} = fullfile(tarDIR, outDIR, strcat(bodyFN,suffixFN,sprintf('_%d.bild',loop)));
    end
    fid = fopen(bildPATH{loop}, 'w');
    
%     % Assign transparency
%     if(isempty(triad))
%         fprintf(fid, '.transparency 0\n');
%     else
%         fprintf(fid, '.transparency 0.5\n');
%     end
    
    for j = 1:size(cylinder{loop},1)
        % Assign colors
        if(isempty(BIN) || cylinder{loop}(j,8)==rTruss)   % without assigned color
            fprintf(fid, '.color gray\n');
        else                                        % with assigned color
            nid = cylinder{loop}(j,1);
            nelem = ceil(j/2);
            if(~isempty(nodeVal))                   % color by node
                binID = BIN(nid,loop);
            elseif(~isempty(elemVal))               % color by element
                binID = BIN(nelem,loop);
            else
                error('Exception.');
            end
            if(binID==0)
                binID = 1;
            end
            fprintf(fid, '.color %f %f %f\n',color(binID,1),color(binID,2),color(binID,3));
        end
        
        % Draw cylinders
        fprintf(fid,'.cylinder %f  %f  %f  %f  %f  %f  %f\n',...
                cylinder{loop}(j,2),cylinder{loop}(j,3),cylinder{loop}(j,4), ...
                cylinder{loop}(j,5),cylinder{loop}(j,6),cylinder{loop}(j,7),cylinder{loop}(j,8));
    end
    
    % Draw alignment elements and related triads
    if(~isempty(triad))
%         fprintf(fid, '.transparency 0\n');
        for j = 1:size(FEModel.alignDSDNA,1)
            nid1 = FEModel.alignDSDNA(j,1);
            nid2 = FEModel.alignDSDNA(j,2);
            triad1 = triad{nid1,loop};
            triad2 = triad{nid2,loop};
            coord1 = pNode(nid1,:,loop);
            coord2 = pNode(nid2,:,loop);
            coord1_a1 = coord1' + triad1(:,1)*lenTriad;
            coord1_a2 = coord1' + triad1(:,2)*lenTriad;
            coord1_a3 = coord1' + triad1(:,3)*lenTriad;
            coord2_a1 = coord2' + triad2(:,1)*lenTriad;
            coord2_a2 = coord2' + triad2(:,2)*lenTriad;
            coord2_a3 = coord2' + triad2(:,3)*lenTriad;
            
            % Draw alignment elements
            if(norm(coord1-coord2)>1e-4)
                fprintf(fid,'.color red\n');
                fprintf(fid,'.cylinder %f  %f  %f  %f  %f  %f  %f\n',...
                        coord1(1), coord1(2), coord1(3), coord2(1), coord2(2), coord2(3), rTruss);
            end
                
            % Draw related triads
            % node 1
            fprintf(fid,'.color firebrick\n');
            fprintf(fid, '.arrow\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n', ...
                    coord1(1), coord1(2), coord1(3), coord1_a1(1), coord1_a1(2), coord1_a1(3), rTriad, rTriad*2);
            fprintf(fid,'.color sea green\n');
            fprintf(fid, '.arrow\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n', ...
                    coord1(1), coord1(2), coord1(3), coord1_a2(1), coord1_a2(2), coord1_a2(3), rTriad, rTriad*2);
            fprintf(fid,'.color steel blue\n');
            fprintf(fid, '.arrow\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n', ...
                    coord1(1), coord1(2), coord1(3), coord1_a3(1), coord1_a3(2), coord1_a3(3), rTriad, rTriad*2);
                
            % node 2
            fprintf(fid,'.color firebrick\n');
            fprintf(fid, '.arrow\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n', ...
                    coord2(1), coord2(2), coord2(3), coord2_a1(1), coord2_a1(2), coord2_a1(3), rTriad, rTriad*2);
            fprintf(fid,'.color sea green\n');
            fprintf(fid, '.arrow\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n', ...
                    coord2(1), coord2(2), coord2(3), coord2_a2(1), coord2_a2(2), coord2_a2(3), rTriad, rTriad*2);
            fprintf(fid,'.color steel blue\n');
            fprintf(fid, '.arrow\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n', ...
                    coord2(1), coord2(2), coord2(3), coord2_a3(1), coord2_a3(2), coord2_a3(3), rTriad, rTriad*2);
        end
    end
    
    fclose(fid);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert the BILD file to three PNG files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pyPATH = fullfile(tarDIR, strcat(bodyFN,suffixFN,'_chimeraScr.py'));
fid = fopen(pyPATH, 'w');
fprintf(fid,'from chimera import runCommand\n');

for loop = 1:nFrame
    if(nFrame==1)
        imageFN1 = fullfile(tarDIR, outDIR, strcat(bodyFN,suffixFN,'_view1.png'));
        imageFN2 = fullfile(tarDIR, outDIR, strcat(bodyFN,suffixFN,'_view2.png'));
        imageFN3 = fullfile(tarDIR, outDIR, strcat(bodyFN,suffixFN,'_view3.png'));
    else
        imageFN1 = fullfile(tarDIR, outDIR, strcat(bodyFN,suffixFN,sprintf('_view1_%03d.png',loop)));
        imageFN2 = fullfile(tarDIR, outDIR, strcat(bodyFN,suffixFN,sprintf('_view2_%03d.png',loop)));
        imageFN3 = fullfile(tarDIR, outDIR, strcat(bodyFN,suffixFN,sprintf('_view3_%03d.png',loop)));
    end
    
    fprintf(fid,'runCommand(''open %s'')\n', strrep(bildPATH{loop},'\','/'));
    fprintf(fid,'runCommand(''windowsize 640 480'')\n');
    fprintf(fid,'runCommand(''preset apply publication 3'')\n');
    fprintf(fid,'runCommand(''window'')\n');
    fprintf(fid,'runCommand(''scale 0.6'')\n');
    
    fprintf(fid,'runCommand(''copy file %s png supersample 3'')\n', strrep(imageFN1,'\','/'));
    fprintf(fid,'runCommand(''wait'')\n');
    fprintf(fid,'runCommand(''turn x 90'')\n');
    fprintf(fid,'runCommand(''copy file %s png supersample 3'')\n', strrep(imageFN2,'\','/'));
    fprintf(fid,'runCommand(''wait'')\n');
    fprintf(fid,'runCommand(''turn x -90'')\n');
    fprintf(fid,'runCommand(''turn y 90'')\n');
    fprintf(fid,'runCommand(''copy file %s png supersample 3'')\n', strrep(imageFN3,'\','/'));
    fprintf(fid,'runCommand(''wait'')\n');
    
    fprintf(fid,'runCommand(''close all'')\n');
end

fprintf(fid,'runCommand(''stop yes'')\n');
fclose(fid);


runChimera = sprintf('%s %s %s',param.chimeraEXE,param.chimeraOPTION,pyPATH);
system(runChimera);

end


function cmap = cmapGen(hue, value, n, option)

if(strcmp(option,'single'))
    assert(numel(hue)==1);
    cmap = hsv2rgb([hue*ones(n,1) linspace(0,1,n)' value*ones(n,1)]);
elseif(strcmp(option,'double'))
    assert(numel(hue)==2);
    assert(mod(n,2)==1);
    m = ceil(n/2);
    cmap1 = hsv2rgb([hue(1)*ones(m,1) linspace(1,0,m)' value*ones(m,1)]);
    cmap1(end,:) = [];
    cmap2 = hsv2rgb([hue(2)*ones(m,1) linspace(0,1,m)' value*ones(m,1)]);
    cmap = [cmap1; cmap2];
else
    error('Unrecognizable option.');
end

end