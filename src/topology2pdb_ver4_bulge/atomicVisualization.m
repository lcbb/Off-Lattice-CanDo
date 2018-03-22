function [] = atomicVisualization(workDIR, jobNAME, visualDIR, strand, sysParam)

jsonPATH = fullfile(workDIR, jobNAME, sprintf('%s.json',jobNAME));
pdbPrefix = fullfile(workDIR, jobNAME);
nmaPrefix = fullfile(workDIR, sprintf('%s_NMA',jobNAME), jobNAME);
movieDIR = fullfile(visualDIR, 'movieAtomic');

pdbPrefix_multimodel = fullfile(visualDIR, sprintf('%s_multimodel',jobNAME));
nmaPrefix_multimodel = fullfile(movieDIR, sprintf('%s_multimodel',jobNAME));

if(exist(visualDIR,'dir')==7)
    rmdir(visualDIR, 's');
end
mkdir(visualDIR);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1. Obtain the RGB value for each strand
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Read the JSON file
% fid = fopen(jsonPATH);
% inString = fscanf(fid,'%c'); 
% fclose(fid);
% data = p_json(inString);
% 
% % Assign colors to each staple strand.
% % First linear strands, then circular strands.
% stapColor = [];
% stapColor.line = [];
% stapColor.loop = [];
% for i = 1:numel(data.vstrands)
%     helixID = data.vstrands{i}.num + 1;     % Which helix?
%     nStap = numel(data.vstrands{i}.stap_colors);
%     for j = 1:nStap
%         stapPos = data.vstrands{i}.stap_colors{j}{1} + 1;   % Which position?
%         currRGB = getRGB(data.vstrands{i}.stap_colors{j}{2});       % Which color?
%         stapID = topology.stapTourMark{helixID}(stapPos,1);
% %         assert(topology.stapTourMark{helixID}(stapPos,2) == 1);
%         if(stapID > 0)
%             stapColor.line(stapID,:) = currRGB;
%         elseif(stapID < 0)
%             stapColor.loop(-stapID,:) = currRGB;
%         else
%             error('Empty strand.');
%         end
%     end
% end
% 
% % Build the color list for all the strands.
% % The scaffold strands always have RGB = [0, 102, 204].
% % See my email to Shawn Douglas titled "caDNAno JSON files"
% % sent on Monday, October 07, 2013 5:50 PM.
% strandColor = repmat([0 102 204], numel(topology.scafTourLine)+numel(topology.scafTourLoop), 1);
% strandColor = cat(1, strandColor, stapColor.line);
% strandColor = cat(1, strandColor, stapColor.loop);

strandColorList = [0 102 204; 184 5 108; 247 67 8; 3 182 162; 247 147 30; 204 0 0; 87 187 0; 0 114 0; 115 0 222];
nColor = size(strandColorList,1);
nStrand = numel(strand);
strandColor = zeros(nStrand,3);
for i = 1:nStrand
    strandColor(i,:) = strandColorList(mod(i-1,nColor)+1,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2. Align the PDB file to the principal axes
%         See CanDo v2.3 script generateRMSFmovie.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(fullfile(workDIR,sprintf('%s_FEModel.mat',jobNAME)), 'FEModel');
load(fullfile(workDIR,sprintf('%s_Displ.mat',jobNAME)), 'displ');
dNode = FEModel.node(:,1:3) + displ(:,1:3);
[pAxes,pNode] = pca(dNode);
if(dot(cross(pAxes(:,1),pAxes(:,2)),pAxes(:,3))<0);
    pAxes(:,3) = -pAxes(:,3);
end
rotAxisAngle = vrrotmat2vec(pAxes');
Rot = [];
Rot.axis = rotAxisAngle(1:3)';
Rot.angle = rotAxisAngle(4);
Rot.center = mean(dNode,1)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3. Make images for the mechanical ground-state structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nModel = pdb2multimodel(sprintf('%s.pdb',pdbPrefix), sprintf('%s.pdb',pdbPrefix_multimodel));
% ssDNA_colorList = generatePDB4ssDNA(deformList, sprintf('%s.pdb',pdbPrefix_multimodel), sprintf('%s_ssDNA.pdb',pdbPrefix_multimodel), [0 102 204], stapColor);
% if(size(strandColor,1)<nModel)
%     strandColor = cat(1,strandColor,repmat([127 127 127],nModel-size(strandColor,1),1));
% end
copyfile(sprintf('%s.pdb',pdbPrefix), sprintf('%s.pdb',pdbPrefix_multimodel));
ssDNA_colorList = generatePDB4ssDNA(sprintf('%s.pdb',pdbPrefix_multimodel), sprintf('%s_ssDNA.pdb',pdbPrefix_multimodel), strand, strandColor);
pdb2png(pdbPrefix_multimodel, pdbPrefix_multimodel, strandColor, ssDNA_colorList, Rot, sysParam);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4. Make a movie for the NMA result
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(sysParam.atomicMovie)
    if(exist(movieDIR,'dir')==7)
        rmdir(movieDIR, 's');
    end
    mkdir(movieDIR);

    % PDB to PNG
    for i = 1:sysParam.nFrame
%         pdb2multimodel(sprintf('%s_%d.pdb',nmaPrefix,i), sprintf('%s_%d.pdb',nmaPrefix_multimodel,i));
%         ssDNA_colorList = generatePDB4ssDNA(deformList, sprintf('%s_%d.pdb',nmaPrefix,i), sprintf('%s_%d_ssDNA.pdb',nmaPrefix_multimodel,i), [0 102 204], stapColor);
        copyfile(sprintf('%s_%d.pdb',nmaPrefix,i), sprintf('%s_%d.pdb',nmaPrefix_multimodel,i));
        ssDNA_colorList = generatePDB4ssDNA(sprintf('%s_%d.pdb',nmaPrefix_multimodel,i), sprintf('%s_%d_ssDNA.pdb',nmaPrefix_multimodel,i), strand, strandColor);
    end
    pdb2png(nmaPrefix_multimodel, nmaPrefix_multimodel, strandColor, ssDNA_colorList, Rot, sysParam, 1:sysParam.nFrame);
    
    % PNG to AVI
    png2avi(nmaPrefix_multimodel, sysParam.nFrame, 10, sysParam);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get RGB value from a single number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RGB = getRGB(int)

RGB = zeros(1,3);
RGB(1) = bitand(bitshift(int,-16),255);
RGB(2) = bitand(bitshift(int,-8),255);
RGB(3) = bitand(int,255);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert the format of a PDB structure: build a model for each chain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function iModel = pdb2multimodel(filenameIn, filenameOut)

fidIn = fopen(filenameIn);
fidOut = fopen(filenameOut,'w');

iModel = 0;
iAtom = 0;
isNewModel = true;

tline = fgetl(fidIn);
while(ischar(tline))
    if(strcmp(tline(1:4),'ATOM'))
        if(isNewModel)
            iModel = iModel + 1;
            fprintf(fidOut, 'MODEL%9d\n', iModel);
            isNewModel = false;
        end
        tline(21:22) = ' A';
%         if(~strcmp(tline(21:22),' A'))
%             tline(21:22) = ' B';
%         end
        iAtom = iAtom + 1;
        tline(6:11) = sprintf('%6d',iAtom);
        fprintf(fidOut, '%s\n', tline);
        %atomSerNo = str2double(tline(6:11));
        %atomSerNo = mod(atomSerNo-1,99999) + 1;
        %fprintf(fidOut, '%s%6.0f%s\n', tline(1:5),atomSerNo,tline(12:end));
    elseif(strcmp(tline(1:3),'TER'))
        tline(21:22) = ' A';
%         if(~strcmp(tline(21:22),' A'))
%             tline(21:22) = ' B';
%         end
        iAtom = iAtom + 1;
        tline(6:11) = sprintf('%6d',iAtom);
        iAtom = 0;
        fprintf(fidOut, '%s\n', tline);
        
        fprintf(fidOut, 'ENDMDL\n');
        isNewModel = true;
    end
    tline = fgetl(fidIn);
end

fclose(fidIn);
fclose(fidOut);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert PDB structures to PNG images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = pdb2png(pdbPrefix, pngPrefix, strandColor, ssDNA_colorList, Rot, sysParam, iFile)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the file names of the PDB and PNG files
% Get the rotation matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(nargin <= 6)
    pdbPATH{1} = sprintf('%s.pdb', pdbPrefix);
    ssDNAPATH{1} = sprintf('%s_ssDNA.pdb', pdbPrefix);
    png1PATH{1} = sprintf('%s_view1.png', pngPrefix);
    png2PATH{1} = sprintf('%s_view2.png', pngPrefix);
    png3PATH{1} = sprintf('%s_view3.png', pngPrefix);
elseif(nargin == 7)
    assert(numel(iFile)>0);
    pdbPATH = cell(numel(iFile),1);
    ssDNAPATH = cell(numel(iFile),1);
    png1PATH = cell(numel(iFile),1);
    png2PATH = cell(numel(iFile),1);
    png3PATH = cell(numel(iFile),1);
    for loop = 1:numel(iFile)
        pdbPATH{loop} = sprintf('%s_%d.pdb', pdbPrefix, loop);
        ssDNAPATH{loop} = sprintf('%s_%d_ssDNA.pdb', pdbPrefix, loop);
        png1PATH{loop} = sprintf('%s_view1_%03d.png', pngPrefix, loop);
        png2PATH{loop} = sprintf('%s_view2_%03d.png', pngPrefix, loop);
        png3PATH{loop} = sprintf('%s_view3_%03d.png', pngPrefix, loop);
    end
else
    error('Incorrect argument format.');
end     

rAxis = Rot.axis;
assert(norm(rAxis)>0);
rAxis = rAxis/norm(rAxis);
rAngle = Rot.angle/pi*180;
rCenter = Rot.center;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the UCSF Chimera script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chimeraScr = sprintf('%s_chimera.py',pngPrefix);
fid = fopen(chimeraScr, 'w');
% Import the Python interface
fprintf(fid,'from chimera import runCommand\n');

for loop = 1:numel(png1PATH)
    % Open the PDB file
    fprintf(fid, 'runCommand(''open %s'')\n', strrep(pdbPATH{loop},'\','/'));
    
    % Set the environment
%     fprintf(fid, 'runCommand(''windowsize 640 480'')\n');
    fprintf(fid, 'runCommand(''windowsize 1280 960'')\n');
    fprintf(fid, 'runCommand(''preset apply publication 3'')\n');
    fprintf(fid, 'runCommand(''window'')\n');
    fprintf(fid, 'runCommand(''scale 0.8'')\n');
    
    % Turn off the original rendering
    fprintf(fid, 'runCommand(''~ribbon'')\n');
    fprintf(fid, 'runCommand(''~display'')\n');
    
    % Use the new rendering
    for i = 1:size(strandColor,1)
        RGB = strandColor(i,:)/255;
        % Check if only one chain
        if (size(strandColor,1) == 1)
            fprintf(fid, 'runCommand(''molmap #0 3'')\n');
            fprintf(fid, 'runCommand(''volume #0 color %f,%f,%f step 1'')\n', RGB(1), RGB(2), RGB(3));
        else
            fprintf(fid, 'runCommand(''molmap #0.%d 3'')\n', i);
            fprintf(fid, 'runCommand(''volume #0.%d color %f,%f,%f step 1'')\n', i, RGB(1), RGB(2), RGB(3));
        end
    end
    
    % Render the ssDNA
    if(~isempty(ssDNA_colorList))
        fprintf(fid, 'runCommand(''open %s'')\n', strrep(ssDNAPATH{loop},'\','/'));
        fprintf(fid, 'runCommand(''~ribbon #1'')\n');
        fprintf(fid, 'runCommand(''~display #1'')\n');
        if(size(ssDNA_colorList,1) == 1)
            fprintf(fid, 'runCommand(''shape tube #1 radius 0.75 color %f,%f,%f'')\n', ssDNA_colorList(i,1)/255, ssDNA_colorList(i,2)/255, ssDNA_colorList(i,3)/255);
        else
            for i = 1:size(ssDNA_colorList,1)
                fprintf(fid, 'runCommand(''shape tube #1.%d radius 0.75 color %f,%f,%f'')\n', i, ssDNA_colorList(i,1)/255, ssDNA_colorList(i,2)/255, ssDNA_colorList(i,3)/255);
            end
        end
    end
    
    % Save as PNG files
    %fprintf(fid, 'runCommand(''turn %f,%f,%f %f center %f,%f,%f'')\n', rAxis(1), rAxis(2), rAxis(3), rAngle, rCenter(1), rCenter(2), rCenter(3));
    fprintf(fid, 'runCommand(''turn %f,%f,%f %f'')\n', rAxis(1), rAxis(2), rAxis(3), rAngle);
    fprintf(fid, 'runCommand(''window'')\n');
    fprintf(fid, 'runCommand(''scale 0.8'')\n');
    fprintf(fid, 'runCommand(''wait'')\n');
    
    fprintf(fid, 'runCommand(''copy file %s png supersample 3'')\n', strrep(png1PATH{loop},'\','/'));
    fprintf(fid, 'runCommand(''wait'')\n');
    %fprintf(fid, 'runCommand(''turn x 90 center %f,%f,%f'')\n', rCenter(1), rCenter(2), rCenter(3));
    fprintf(fid, 'runCommand(''turn x 90'')\n');
    fprintf(fid, 'runCommand(''copy file %s png supersample 3'')\n', strrep(png2PATH{loop},'\','/'));
    fprintf(fid, 'runCommand(''wait'')\n');
    %fprintf(fid, 'runCommand(''turn x -90 center %f,%f,%f'')\n', rCenter(1), rCenter(2), rCenter(3));
    fprintf(fid, 'runCommand(''turn x -90'')\n');
    %fprintf(fid, 'runCommand(''turn y 90 center %f,%f,%f'')\n', rCenter(1), rCenter(2), rCenter(3));
    fprintf(fid, 'runCommand(''turn y 90'')\n');
    fprintf(fid, 'runCommand(''copy file %s png supersample 3'')\n', strrep(png3PATH{loop},'\','/'));
    fprintf(fid, 'runCommand(''wait'')\n');
    %fprintf(fid, 'runCommand(''turn y -90 center %f,%f,%f'')\n', rCenter(1), rCenter(2), rCenter(3));
    fprintf(fid, 'runCommand(''turn y -90'')\n');
    fprintf(fid, 'runCommand(''close all'')\n');
end
fprintf(fid, 'runCommand(''stop yes'')\n');
fclose(fid);

runChimera = sprintf('%s %s %s',sysParam.chimeraEXE, sysParam.chimeraOPTION, chimeraScr);
system(runChimera);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert PNG images to AVI movies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function png2avi(pngPATHpart, nSnapshot, framerate, sysParam)

% Get the directory and file names for the movies
[movieDIR,bodyFN] = fileparts(strcat(pngPATHpart,'.pdb'));
% DEBUG
fprintf('DEBUG MADE IT TO png2avi.m function in atomicVisualization');
disp(movieDIR);
disp(sysParam.ffmpegBIN);

% Generate VideoMach scripts to convert images to movies
%project1 = sprintf('%s_view1.vmp', pngPATHpart);
%project2 = sprintf('%s_view2.vmp', pngPATHpart);
%project3 = sprintf('%s_view3.vmp', pngPATHpart);
%fid1 = fopen(project1,'w');
%fid2 = fopen(project2,'w');
%fid3 = fopen(project3,'w');
%writeVideoMachProject(fid1,movieDIR,sprintf('%s_view1',bodyFN),nSnapshot,'avi',framerate,sysParam);
%writeVideoMachProject(fid2,movieDIR,sprintf('%s_view2',bodyFN),nSnapshot,'avi',framerate,sysParam);
%writeVideoMachProject(fid3,movieDIR,sprintf('%s_view3',bodyFN),nSnapshot,'avi',framerate,sysParam);
%fclose(fid1);
%fclose(fid2);
%fclose(fid3);

% No Videomach on Linux
% Run VideoMach
%runVideoMach1 = sprintf('%s /openproject="%s" /start /exit', sysParam.videomachEXE, project1);
%runVideoMach2 = sprintf('%s /openproject="%s" /start /exit', sysParam.videomachEXE, project2);
%runVideoMach3 = sprintf('%s /openproject="%s" /start /exit', sysParam.videomachEXE, project3);
%system(runVideoMach1);
%system(runVideoMach2);
%system(runVideoMach3);

% FFMPEG
cd(movieDIR);
% Run ffmpeg
cd(movieDIR);
% ffmpeg strings
ffmCmd1 = sprintf('-thread_queue_size 512 -framerate 10 -i multimodel_view1_%%03d.png -i /scratch/cando/CanDo/Onlattice/CanDo-logo-new-small.png -c:v libx264 -filter_complex "overlay=230:main_h-overlay_h-10: x=230:y=410" video1.mp4');
ffmCmd2 = sprintf('-thread_queue_size 512 -framerate 10 -i multimodel_view2_%%03d.png -i /scratch/cando/CanDo/Onlattice/CanDo-logo-new-small.png -c:v libx264 -filter_complex "overlay=230:main_h-overlay_h-10: x=230:y=410" video2.mp4');
ffmCmd3 = sprintf('-thread_queue_size 512 -framerate 10 -i multimodel_view3_%%03d.png -i /scratch/cando/CanDo/Onlattice/CanDo-logo-new-small.png -c:v libx264 -filter_complex "overlay=230:main_h-overlay_h-10: x=230:y=410" video3.mp4');
runVideoFFM1 = sprintf('%s %s', sysParam.ffmpegBIN,ffmCmd1);
runVideoFFM2 = sprintf('%s %s', sysParam.ffmpegBIN,ffmCmd2);
runVideoFFM3 = sprintf('%s %s', sysParam.ffmpegBIN,ffmCmd3);
system(runVideoFFM1);
system(runVideoFFM2);
system(runVideoFFM3);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate a VideoMach project file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = writeVideoMachProject(fid,movieDIR,bodyFN,nSnapshot,outputFormat,framerate,param)
%framerate = param.framerate;
%colorbarFN = param.colorbarFN;
logoFN = param.logoFN;
% writing header
fprintf(fid,'<Header>\n');
fprintf(fid,'	<BlockVersion type="int">314</>\n');
fprintf(fid,'	<Id type="string">VideoMach#20Project</>\n');
fprintf(fid,'	<Creator type="string">VideoMach#205.8.3</>\n');
fprintf(fid,'</Header>\n');
% writing file list
fprintf(fid,'<FileList>\n');
fprintf(fid,'	<Video>\n');
fprintf(fid,'		<BlockVersion type="int">314</>\n');
for i=0:nSnapshot
    fprintf(fid,'		<Item%d>\n',i+1);
    fprintf(fid,'			<UniqueId type="int">%d</>\n',i);
    fprintf(fid,'			<FileName type="string">%s_%03d.png</>\n',fullfile(movieDIR,bodyFN),i);
    fprintf(fid,'			<StartQuant type="int">0</>\n');
    fprintf(fid,'			<EndQuant type="int">0</>\n');
    fprintf(fid,'			<MaxQuants type="int">1</>\n');
    fprintf(fid,'			<FrameRate type="float">%d</>\n',framerate);
    fprintf(fid,'			<ReadFwd type="bool">T</>\n');
    fprintf(fid,'			<LinkId type="int">2147483647</>\n');
    fprintf(fid,'		</Item%d>\n',i+1);
end
fprintf(fid,'	</Video>\n');
fprintf(fid,'</FileList>\n');
fprintf(fid,'<OpticFilters>\n');
fprintf(fid,'	<BlockVersion type="int">314</>\n');
% fprintf(fid,'	<Item1>\n');
% fprintf(fid,'		<Id type="string">Overlay</>\n');
% fprintf(fid,'		<Data1>\n');
% fprintf(fid,'			<FileName type="string">%s</>\n',colorbarFN);
% fprintf(fid,'			<SeqMode type="int">1024</>\n');
% fprintf(fid,'			<Corner type="int">1634497125</>\n');
% fprintf(fid,'			<OffsX type="int">400</>\n');
% fprintf(fid,'			<OffsY type="int">400</>\n');
% fprintf(fid,'			<LoopMode type="int">1</>\n');
% fprintf(fid,'			<LoopCount type="int">0</>\n');
% fprintf(fid,'			<TranspMode type="int">2</>\n');
% fprintf(fid,'			<TranspColor type="int">0</>\n');
% fprintf(fid,'			<AntiAliasing type="bool">T</>\n');
% fprintf(fid,'		</Data1>\n');
% fprintf(fid,'	</Item1>\n');
% fprintf(fid,'	<Item2>\n');
% fprintf(fid,'		<Id type="string">AddText</>\n');
% fprintf(fid,'		<Data1>\n');
% fprintf(fid,'			<FontName type="string">Arial</>\n');
% fprintf(fid,'			<FontHeight type="int">25</>\n');
% fprintf(fid,'			<CharSet type="word">1</>\n');
% fprintf(fid,'			<TextColor type="int">16777215</>\n');
% fprintf(fid,'			<FontBold type="bool">F</>\n');
% fprintf(fid,'			<FontItalic type="bool">F</>\n');
% fprintf(fid,'			<FontUnderline type="bool">F</>\n');
% fprintf(fid,'			<FontStrikeout type="bool">F</>\n');
% fprintf(fid,'			<CharAngle type="float">0</>\n');
% fprintf(fid,'			<TextAngle type="float">0</>\n');
% fprintf(fid,'			<ShadowEnabled type="bool">T</>\n');
% fprintf(fid,'			<ShadowColor type="int">0</>\n');
% fprintf(fid,'			<ShadowBlur type="int">2</>\n');
% fprintf(fid,'			<ShadowOpacity type="int">100</>\n');
% fprintf(fid,'			<ShadowOffsetX type="int">2</>\n');
% fprintf(fid,'			<ShadowOffsetY type="int">2</>\n');
% fprintf(fid,'			<TextString type="string">%1.1f</>\n',minRMSFnm);
% fprintf(fid,'			<XPos type="int">360</>\n');
% fprintf(fid,'			<YPos type="int">445</>\n');
% fprintf(fid,'			<TimeOffset type="float">0</>\n');
% fprintf(fid,'			<TimeScale type="float">1</>\n');
% fprintf(fid,'			<GUISameFont type="bool">T</>\n');
% fprintf(fid,'		</Data1>\n');
% fprintf(fid,'	</Item2>\n');
% fprintf(fid,'	<Item3>\n');
% fprintf(fid,'		<Id type="string">AddText</>\n');
% fprintf(fid,'		<Data1>\n');
% fprintf(fid,'			<FontName type="string">Arial</>\n');
% fprintf(fid,'			<FontHeight type="int">25</>\n');
% fprintf(fid,'			<CharSet type="word">1</>\n');
% fprintf(fid,'			<TextColor type="int">16777215</>\n');
% fprintf(fid,'			<FontBold type="bool">F</>\n');
% fprintf(fid,'			<FontItalic type="bool">F</>\n');
% fprintf(fid,'			<FontUnderline type="bool">F</>\n');
% fprintf(fid,'			<FontStrikeout type="bool">F</>\n');
% fprintf(fid,'			<CharAngle type="float">0</>\n');
% fprintf(fid,'			<TextAngle type="float">0</>\n');
% fprintf(fid,'			<ShadowEnabled type="bool">T</>\n');
% fprintf(fid,'			<ShadowColor type="int">0</>\n');
% fprintf(fid,'			<ShadowBlur type="int">2</>\n');
% fprintf(fid,'			<ShadowOpacity type="int">100</>\n');
% fprintf(fid,'			<ShadowOffsetX type="int">2</>\n');
% fprintf(fid,'			<ShadowOffsetY type="int">2</>\n');
% fprintf(fid,'			<TextString type="string">#3E%1.1f#20nm</>\n',maxRMSFnm);
% fprintf(fid,'			<XPos type="int">560</>\n');
% fprintf(fid,'			<YPos type="int">445</>\n');
% fprintf(fid,'			<TimeOffset type="float">0</>\n');
% fprintf(fid,'			<TimeScale type="float">1</>\n');
% fprintf(fid,'			<GUISameFont type="bool">T</>\n');
% fprintf(fid,'		</Data1>\n');
% fprintf(fid,'	</Item3>\n');
fprintf(fid,'	<Item4>\n');
fprintf(fid,'		<Id type="string">Overlay</>\n');
fprintf(fid,'		<Data1>\n');
fprintf(fid,'			<FileName type="string">%s</>\n',logoFN);
fprintf(fid,'			<SeqMode type="int">1024</>\n');
fprintf(fid,'			<Corner type="int">88</>\n');
fprintf(fid,'			<OffsX type="int">0</>\n');
fprintf(fid,'			<OffsY type="int">400</>\n');
fprintf(fid,'			<LoopMode type="int">1</>\n');
fprintf(fid,'			<LoopCount type="int">0</>\n');
fprintf(fid,'			<TranspMode type="int">0</>\n');
fprintf(fid,'			<TranspColor type="int">15395562</>\n');
fprintf(fid,'			<AntiAliasing type="bool">T</>\n');
fprintf(fid,'		</Data1>\n');
fprintf(fid,'	</Item4>\n');
fprintf(fid,'	<Item5>\n');
fprintf(fid,'		<Id type="string">AddText</>\n');
fprintf(fid,'		<Data1>\n');
fprintf(fid,'			<FontName type="string">Arial</>\n');
fprintf(fid,'			<FontHeight type="int">22</>\n');
fprintf(fid,'			<CharSet type="word">1</>\n');
fprintf(fid,'			<TextColor type="int">16777215</>\n');
fprintf(fid,'			<FontBold type="bool">F</>\n');
fprintf(fid,'			<FontItalic type="bool">F</>\n');
fprintf(fid,'			<FontUnderline type="bool">F</>\n');
fprintf(fid,'			<FontStrikeout type="bool">F</>\n');
fprintf(fid,'			<CharAngle type="float">0</>\n');
fprintf(fid,'			<TextAngle type="float">0</>\n');
fprintf(fid,'			<ShadowEnabled type="bool">T</>\n');
fprintf(fid,'			<ShadowColor type="int">0</>\n');
fprintf(fid,'			<ShadowBlur type="int">1</>\n');
fprintf(fid,'			<ShadowOpacity type="int">100</>\n');
fprintf(fid,'			<ShadowOffsetX type="int">2</>\n');
fprintf(fid,'			<ShadowOffsetY type="int">2</>\n');
fprintf(fid,'			<TextString type="string">http://cando-dna-origami.org</>\n');
fprintf(fid,'			<XPos type="int">50</>\n');
fprintf(fid,'			<YPos type="int">450</>\n');
fprintf(fid,'			<TimeOffset type="float">0</>\n');
fprintf(fid,'			<TimeScale type="float">1</>\n');
fprintf(fid,'			<GUISameFont type="bool">T</>\n');
fprintf(fid,'		</Data1>\n');
fprintf(fid,'	</Item5>\n');
fprintf(fid,'</OpticFilters>\n');
% writing output format
if strcmp(outputFormat,'avi')
    fprintf(fid,'<OutputSettings>\n');
    fprintf(fid,'	<Common>\n');   
    fprintf(fid,'		<BlockVersion type="int">314</>\n');
    fprintf(fid,'		<OutputMode type="int">0</>\n');
    fprintf(fid,'		<SplitSize type="float">2147483648</>\n');
    fprintf(fid,'		<SaveInfo>\n');
    fprintf(fid,'			<Version type="word">2</>\n');
    fprintf(fid,'			<ffComments>\n');
    fprintf(fid,'			</ffComments>\n');
    fprintf(fid,'		</SaveInfo>\n');
    fprintf(fid,'	</Common>\n');
    fprintf(fid,'	<Video>\n');
    fprintf(fid,'		<BlockVersion type="int">500</>\n');
    fprintf(fid,'		<FileName type="string">%s.avi</>\n',fullfile(movieDIR,strcat('fluctuation',bodyFN)));
    fprintf(fid,'		<FrameRateCustom type="bool">F</>\n');
    fprintf(fid,'		<FrameRate type="float">10</>\n');
    fprintf(fid,'		<KeepDuration type="bool">T</>\n');
    fprintf(fid,'		<ColorDepthCustom type="bool">F</>\n');
    fprintf(fid,'		<ColorDepth type="int">24</>\n');
    fprintf(fid,'		<SinglePalette type="bool">T</>\n');
    fprintf(fid,'		<Dithering type="int">2</>\n');
    fprintf(fid,'		<GrayScale type="bool">F</>\n');
    fprintf(fid,'		<ResizeCustom type="bool">F</>\n');
    fprintf(fid,'		<ResizeFolder>\n');
    fprintf(fid,'			<ResX type="int">0</>\n');
    fprintf(fid,'			<ResY type="int">0</>\n');
    fprintf(fid,'			<Bilinear type="bool">T</>\n');
    fprintf(fid,'			<EnlargeMode type="int">0</>\n');
    fprintf(fid,'			<ShrinkMode type="int">0</>\n');
    fprintf(fid,'			<BorderColor type="int">0</>\n');
    fprintf(fid,'		</ResizeFolder>\n');
    fprintf(fid,'		<PrintResCustom type="bool">F</>\n');
    fprintf(fid,'		<PrintResX type="int">0</>\n');
    fprintf(fid,'		<PrintResY type="int">0</>\n');
    fprintf(fid,'		<VideoFormatInfo type="string">Cinepak#20Codec#20by#20Radius</>\n');
    fprintf(fid,'		<FileFormatIndex type="int">5</>\n');
    fprintf(fid,'		<FileFormatFolder>\n');
    fprintf(fid,'			<vFileFmtId type="int">6</>\n');
    fprintf(fid,'			<vCodecId type="int">1684633187</>\n');
    fprintf(fid,'			<vCodecDataSize type="int">0</>\n');
    fprintf(fid,'			<vFPS type="float">10</>\n');
    fprintf(fid,'			<vKeepDuration type="bool">T</>\n');
    fprintf(fid,'			<vBitRate type="int">0</>\n');
    fprintf(fid,'			<vBitRateMax type="int">0</>\n');
    fprintf(fid,'			<vBitRateMin type="int">0</>\n');
    fprintf(fid,'			<vKeyFrameRate type="int">0</>\n');
    fprintf(fid,'			<vQuality type="float">100</>\n');
    fprintf(fid,'			<vAsr type="int">0</>\n');
    fprintf(fid,'			<vDpiX type="int">0</>\n');
    fprintf(fid,'			<vDpiY type="int">0</>\n');
    fprintf(fid,'			<vOutResX type="int">640</>\n');
    fprintf(fid,'			<vOutResY type="int">480</>\n');
    fprintf(fid,'			<vOutBpp type="int">24</>\n');
    fprintf(fid,'			<vCodecName type="string">Cinepak#20Codec#20by#20Radius</>\n');
    fprintf(fid,'		</FileFormatFolder>\n');
    fprintf(fid,'	</Video>\n');
    fprintf(fid,'	<Audio>\n');
    fprintf(fid,'		<BlockVersion type="int">500</>\n');
    fprintf(fid,'		<FileName type="string">%s.avi</>\n',fullfile(movieDIR,bodyFN));
    fprintf(fid,'		<SampleRateCustom type="bool">F</>\n');
    fprintf(fid,'		<SampleRate type="word">44100</>\n');
    fprintf(fid,'		<ChannelsCustom type="bool">F</>\n');
    fprintf(fid,'		<Channels type="word">2</>\n');
    fprintf(fid,'		<SampleResCustom type="bool">F</>\n');
    fprintf(fid,'		<SampleRes type="word">16</>\n');
    fprintf(fid,'		<AudioFormatInfo type="string">PCM,#20uncompressed</>\n');
    fprintf(fid,'		<FileFormatIndex type="int">0</>\n');
    fprintf(fid,'		<FileFormatFolder>\n');
    fprintf(fid,'			<aFileFmtId type="int">6</>\n');
    fprintf(fid,'			<aCodecId type="word">1</>\n');
    fprintf(fid,'			<aWaveFmtOut type="bin" size="18">0100020044AC000010B10200040010000000</>\n');
    fprintf(fid,'			<aBitRate type="word">0</>\n');
    fprintf(fid,'			<aQuality type="float">0</>\n');
    fprintf(fid,'			<aCodecName type="string">PCM,#20uncompressed</>\n');
    fprintf(fid,'		</FileFormatFolder>\n');
    fprintf(fid,'	</Audio>\n');
    fprintf(fid,'</OutputSettings>    \n');
elseif strcmp(outputFormat,'png')
    fprintf(fid,'<OutputSettings>n');
    fprintf(fid,'	<Common>n');
    fprintf(fid,'		<BlockVersion type="int">314</>n');
    fprintf(fid,'		<OutputMode type="int">0</>n');
    fprintf(fid,'		<SplitSize type="float">1099511627776000000</>n');
    fprintf(fid,'		<SaveInfo>n');
    fprintf(fid,'			<Version type="word">2</>n');
    fprintf(fid,'			<ffComments>n');
    fprintf(fid,'			</ffComments>n');
    fprintf(fid,'		</SaveInfo>n');
    fprintf(fid,'	</Common>n');
    fprintf(fid,'	<Video>n');
    fprintf(fid,'		<BlockVersion type="int">500</>n');
    fprintf(fid,'		<FileName type="string">%s_averaged.png</>n',fullfile(movieDIR,bodyFN));
    fprintf(fid,'		<FrameRateCustom type="bool">T</>n');
    fprintf(fid,'		<FrameRate type="float">10</>n');
    fprintf(fid,'		<KeepDuration type="bool">F</>n');
    fprintf(fid,'		<ColorDepthCustom type="bool">F</>n');
    fprintf(fid,'		<ColorDepth type="int">24</>n');
    fprintf(fid,'		<SinglePalette type="bool">T</>n');
    fprintf(fid,'		<Dithering type="int">2</>n');
    fprintf(fid,'		<GrayScale type="bool">F</>n');
    fprintf(fid,'		<ResizeCustom type="bool">F</>n');
    fprintf(fid,'		<ResizeFolder>n');
    fprintf(fid,'			<ResX type="int">0</>n');
    fprintf(fid,'			<ResY type="int">0</>n');
    fprintf(fid,'			<Bilinear type="bool">T</>n');
    fprintf(fid,'			<EnlargeMode type="int">0</>n');
    fprintf(fid,'			<ShrinkMode type="int">0</>n');
    fprintf(fid,'			<BorderColor type="int">0</>n');
    fprintf(fid,'		</ResizeFolder>n');
    fprintf(fid,'		<PrintResCustom type="bool">F</>n');
    fprintf(fid,'		<PrintResX type="int">0</>n');
    fprintf(fid,'		<PrintResY type="int">0</>n');
    fprintf(fid,'		<VideoFormatInfo type="string">PNG#20Compression#20Level:#207#20of#209</>n');
    fprintf(fid,'		<FileFormatIndex type="int">8</>n');
    fprintf(fid,'		<FileFormatFolder>n');
    fprintf(fid,'			<vFileFmtId type="int">9</>n');
    fprintf(fid,'			<vCodecId type="int">7</>n');
    fprintf(fid,'			<vCodecDataSize type="int">0</>n');
    fprintf(fid,'			<vFPS type="float">10</>n');
    fprintf(fid,'			<vKeepDuration type="bool">F</>n');
    fprintf(fid,'			<vBitRate type="int">0</>n');
    fprintf(fid,'			<vBitRateMax type="int">0</>n');
    fprintf(fid,'			<vBitRateMin type="int">0</>n');
    fprintf(fid,'			<vKeyFrameRate type="int">0</>n');
    fprintf(fid,'			<vQuality type="float">0</>n');
    fprintf(fid,'			<vAsr type="int">0</>n');
    fprintf(fid,'			<vDpiX type="int">0</>n');
    fprintf(fid,'			<vDpiY type="int">0</>n');
    fprintf(fid,'			<vOutResX type="int">640</>n');
    fprintf(fid,'			<vOutResY type="int">480</>n');
    fprintf(fid,'			<vOutBpp type="int">24</>n');
    fprintf(fid,'			<vCodecName type="string">PNG#20Compression#20Level:#207#20of#209</>n');
    fprintf(fid,'		</FileFormatFolder>n');
    fprintf(fid,'	</Video>n');
    fprintf(fid,'	<Audio>n');
    fprintf(fid,'		<BlockVersion type="int">500</>n');
    fprintf(fid,'		<FileName type="string">%s_averaged.png</>n',fullfile(movieDIR,bodyFN));
    fprintf(fid,'		<SampleRateCustom type="bool">F</>n');
    fprintf(fid,'		<SampleRate type="word">44100</>n');
    fprintf(fid,'		<ChannelsCustom type="bool">F</>n');
    fprintf(fid,'		<Channels type="word">2</>n');
    fprintf(fid,'		<SampleResCustom type="bool">F</>n');
    fprintf(fid,'		<SampleRes type="word">16</>n');
    fprintf(fid,'	</Audio>n');
    fprintf(fid,'</OutputSettings>\n');
elseif strcmp(outputFormat,'gif')
    fprintf(fid,'<OutputSettings>\n');
    fprintf(fid,'	<Common>\n');   
    fprintf(fid,'		<BlockVersion type="int">314</>\n');
    fprintf(fid,'		<OutputMode type="int">0</>\n');
    fprintf(fid,'		<SplitSize type="float">2147483648</>\n');
    fprintf(fid,'		<SaveInfo>\n');
    fprintf(fid,'			<Version type="word">2</>\n');
    fprintf(fid,'			<ffComments>\n');
    fprintf(fid,'			</ffComments>\n');
    fprintf(fid,'		</SaveInfo>\n');
    fprintf(fid,'	</Common>\n');
    fprintf(fid,'	<Video>\n');
    fprintf(fid,'		<BlockVersion type="int">500</>\n');
    fprintf(fid,'		<FileName type="string">%s.gif</>\n',fullfile(movieDIR,strcat('fluctuation',bodyFN)));
    fprintf(fid,'		<FrameRateCustom type="bool">F</>\n');
    fprintf(fid,'		<FrameRate type="float">10</>\n');
    fprintf(fid,'		<KeepDuration type="bool">T</>\n');
    fprintf(fid,'		<ColorDepthCustom type="bool">F</>\n');
    fprintf(fid,'		<ColorDepth type="int">24</>\n');
    fprintf(fid,'		<SinglePalette type="bool">T</>\n');
    fprintf(fid,'		<Dithering type="int">2</>\n');
    fprintf(fid,'		<GrayScale type="bool">F</>\n');
    fprintf(fid,'		<ResizeCustom type="bool">F</>\n');
    fprintf(fid,'		<ResizeFolder>\n');
    fprintf(fid,'			<ResX type="int">0</>\n');
    fprintf(fid,'			<ResY type="int">0</>\n');
    fprintf(fid,'			<Bilinear type="bool">T</>\n');
    fprintf(fid,'			<EnlargeMode type="int">0</>\n');
    fprintf(fid,'			<ShrinkMode type="int">0</>\n');
    fprintf(fid,'			<BorderColor type="int">0</>\n');
    fprintf(fid,'		</ResizeFolder>\n');
    fprintf(fid,'		<PrintResCustom type="bool">F</>\n');
    fprintf(fid,'		<PrintResX type="int">0</>\n');
    fprintf(fid,'		<PrintResY type="int">0</>\n');
    fprintf(fid,'		<VideoFormatInfo type="string">Cinepak#20Codec#20by#20Radius</>\n');
    fprintf(fid,'		<FileFormatIndex type="int">5</>\n');
    fprintf(fid,'		<FileFormatFolder>\n');
    fprintf(fid,'			<vFileFmtId type="int">6</>\n');
    fprintf(fid,'			<vCodecId type="int">1684633187</>\n');
    fprintf(fid,'			<vCodecDataSize type="int">0</>\n');
    fprintf(fid,'			<vFPS type="float">10</>\n');
    fprintf(fid,'			<vKeepDuration type="bool">T</>\n');
    fprintf(fid,'			<vBitRate type="int">0</>\n');
    fprintf(fid,'			<vBitRateMax type="int">0</>\n');
    fprintf(fid,'			<vBitRateMin type="int">0</>\n');
    fprintf(fid,'			<vKeyFrameRate type="int">0</>\n');
    fprintf(fid,'			<vQuality type="float">100</>\n');
    fprintf(fid,'			<vAsr type="int">0</>\n');
    fprintf(fid,'			<vDpiX type="int">0</>\n');
    fprintf(fid,'			<vDpiY type="int">0</>\n');
    fprintf(fid,'			<vOutResX type="int">640</>\n');
    fprintf(fid,'			<vOutResY type="int">480</>\n');
    fprintf(fid,'			<vOutBpp type="int">24</>\n');
    fprintf(fid,'			<vCodecName type="string">Cinepak#20Codec#20by#20Radius</>\n');
    fprintf(fid,'		</FileFormatFolder>\n');
    fprintf(fid,'	</Video>\n');
    fprintf(fid,'	<Audio>\n');
    fprintf(fid,'		<BlockVersion type="int">500</>\n');
    fprintf(fid,'		<FileName type="string">%s.gif</>\n',fullfile(movieDIR,strcat('fluctuation',bodyFN)));
    fprintf(fid,'		<SampleRateCustom type="bool">F</>\n');
    fprintf(fid,'		<SampleRate type="word">44100</>\n');
    fprintf(fid,'		<ChannelsCustom type="bool">F</>\n');
    fprintf(fid,'		<Channels type="word">2</>\n');
    fprintf(fid,'		<SampleResCustom type="bool">F</>\n');
    fprintf(fid,'		<SampleRes type="word">16</>\n');
    fprintf(fid,'		<AudioFormatInfo type="string">PCM,#20uncompressed</>\n');
    fprintf(fid,'		<FileFormatIndex type="int">0</>\n');
    fprintf(fid,'		<FileFormatFolder>\n');
    fprintf(fid,'			<aFileFmtId type="int">6</>\n');
    fprintf(fid,'			<aCodecId type="word">1</>\n');
    fprintf(fid,'			<aWaveFmtOut type="bin" size="18">0100020044AC000010B10200040010000000</>\n');
    fprintf(fid,'			<aBitRate type="word">0</>\n');
    fprintf(fid,'			<aQuality type="float">0</>\n');
    fprintf(fid,'			<aCodecName type="string">PCM,#20uncompressed</>\n');
    fprintf(fid,'		</FileFormatFolder>\n');
    fprintf(fid,'	</Audio>\n');
    fprintf(fid,'</OutputSettings>    \n');
else
    error('Undefined output format for VideoMach!!!\n');
end

end
