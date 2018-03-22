function pdbFinal = pdbGenerate(strand)

%% Parameters
refName_A = 'AAA.pdb';
refName_G = 'GGG.pdb';
refName_C = 'CCC.pdb';
refName_T = 'TTT.pdb';


%% Read 8 reference structures:
% - A running forward (32 atoms)
% - G running forward (33 atoms)
% - C running forward (30 atoms)
% - T running forward (32 atoms)
% - T running backward
% - C running backward
% - G running backward
% - A running backward
[ref.A_for,ref.T_rev] = readReference(refName_A);
[ref.G_for,ref.C_rev] = readReference(refName_G);
[ref.C_for,ref.G_rev] = readReference(refName_C);
[ref.T_for,ref.A_rev] = readReference(refName_T);

nStrand = numel(strand);

pdbFinal = [];
pdbFinal.Model.Atom = [];
pdbFinal.Model.Terminal = [];

firstAtomSerNo = 1;     % serial number for the first atom in a tour
for i = 1:nStrand
    [pdbStruct, firstAtomSerNo] = pdbTour(strand(i), ref, firstAtomSerNo, 0);
    pdbFinal.Model.Atom = [pdbFinal.Model.Atom, pdbStruct.Model.Atom];
    pdbFinal.Model.Terminal = [pdbFinal.Model.Terminal, pdbStruct.Model.Terminal];
end




%% Rotate and translate reference structures:
% nScafTourLine = length(deformList.scafTourLine);
% nScafTourLoop = length(deformList.scafTourLoop);
% nStapTourLine = length(deformList.stapTourLine);
% nStapTourLoop = length(deformList.stapTourLoop);
% 
% % PDB structures
% pdbScafTourLine = cell(nScafTourLine,1);
% pdbScafTourLoop = cell(nScafTourLoop,1);
% pdbStapTourLine = cell(nStapTourLine,1);
% pdbStapTourLoop = cell(nStapTourLoop,1);
% 
% pdbFinal = [];
% pdbFinal.Model.Atom = [];
% pdbFinal.Model.Terminal = [];
% 
% % Open-loop scaffold
% firstAtomSerNo = 1;     % serial number for the first atom in a tour
% for i = 1:nScafTourLine
%     [pdbScafTourLine{i}, firstAtomSerNo] = pdbTour(deformList.scafTourLine{i}, ...
%                                                    deformList.scafTourLineR{i}, ...
%                                                    deformList.scafTourLineD{i}, ...
%                                                    deformList.scafTourLineSeq{i}, ...
%                                                    deformList.scafTourLineIndex{i}, ...
%                                                    ref, firstAtomSerNo, 0);
%     pdbFinal.Model.Atom = [pdbFinal.Model.Atom, pdbScafTourLine{i}.Model.Atom];
%     pdbFinal.Model.Terminal = [pdbFinal.Model.Terminal, pdbScafTourLine{i}.Model.Terminal];
% end
% % Close-loop scaffold
% for i = 1:nScafTourLoop
%     [pdbScafTourLoop{i}, firstAtomSerNo] = pdbTour(deformList.scafTourLoop{i}, ...
%                                                    deformList.scafTourLoopR{i}, ...
%                                                    deformList.scafTourLoopD{i}, ...
%                                                    deformList.scafTourLoopSeq{i}, ...
%                                                    deformList.scafTourLoopIndex{i}, ...
%                                                    ref, firstAtomSerNo, 0);
%     pdbFinal.Model.Atom = [pdbFinal.Model.Atom, pdbScafTourLoop{i}.Model.Atom];
%     pdbFinal.Model.Terminal = [pdbFinal.Model.Terminal, pdbScafTourLoop{i}.Model.Terminal];
% end
% % Open-loop staple
% for i = 1:nStapTourLine
%     [pdbStapTourLine{i}, firstAtomSerNo] = pdbTour(deformList.stapTourLine{i}, ...
%                                                    deformList.stapTourLineR{i}, ...
%                                                    deformList.stapTourLineD{i}, ...
%                                                    deformList.stapTourLineSeq{i}, ...
%                                                    deformList.stapTourLineIndex{i}, ...
%                                                    ref, firstAtomSerNo, 1);
%     pdbFinal.Model.Atom = [pdbFinal.Model.Atom, pdbStapTourLine{i}.Model.Atom];
%     pdbFinal.Model.Terminal = [pdbFinal.Model.Terminal, pdbStapTourLine{i}.Model.Terminal];
% end
% % Close-loop staple
% for i = 1:nStapTourLoop
%     [pdbStapTourLoop{i}, firstAtomSerNo] = pdbTour(deformList.stapTourLoop{i}, ...
%                                                    deformList.stapTourLoopR{i}, ...
%                                                    deformList.stapTourLoopD{i}, ...
%                                                    deformList.stapTourLoopSeq{i}, ...
%                                                    deformList.stapTourLoopIndex{i}, ...
%                                                    ref, firstAtomSerNo, 1);
%     pdbFinal.Model.Atom = [pdbFinal.Model.Atom, pdbStapTourLoop{i}.Model.Atom];
%     pdbFinal.Model.Terminal = [pdbFinal.Model.Terminal, pdbStapTourLoop{i}.Model.Terminal];
% end

end


% Build the PDB structure for one tour. Flag = 0--scaffold or 1--staple
function [pdbStruct, firstAtomSerNo_out] = pdbTour(strand, ...
                                                   refStruct, firstAtomSerNo_in, flag)
assert(numel(strand) == 1);
                                               
firstAtomSerNo_out = firstAtomSerNo_in;
pdbStruct = [];
maxAtomPerBase = 40;

nBase = numel(strand.tour);
pdbStruct.Model.Atom = repmat(refStruct.A_for.Model.Atom(1),1,maxAtomPerBase*nBase);

% Rotate about Y_axis about 180 deg
Ry180 = Ry(180);

% Build the PDB structure for each base in the tour
atomIndex = 1;
for i = 1:nBase
    assert((isempty(strand.R{i})&&isempty(strand.d{i})) || (~isempty(strand.R{i})&&~isempty(strand.d{i})));
    if(~isempty(strand.R{i}))
        % Which base?
        b = strand.seq(i);
        % Rotation matrix
        R = strand.R{i}*Ry180;
        %     if(strand.isMain(i))
        %         R = R*Ry180;
        %     end
        % Translation
        D = strand.d{i};
        
        % Read a reference base
        if(strand.isMain(i) == true)        % scaffold strand
            if('A'==b || 'a'==b)
                currentStruct = refStruct.A_for;
            elseif('G'==b || 'g'==b)
                currentStruct = refStruct.G_for;
            elseif('C'==b || 'c'==b)
                currentStruct = refStruct.C_for;
            elseif('T'==b || 't'==b)
                currentStruct = refStruct.T_for;
            else
                error('Illegal base');
            end
        elseif(strand.isMain(i) == false)   % staple strand
            if('A'==b || 'a'==b)
                currentStruct = refStruct.A_rev;
            elseif('G'==b || 'g'==b)
                currentStruct = refStruct.G_rev;
            elseif('C'==b || 'c'==b)
                currentStruct = refStruct.C_rev;
            elseif('T'==b || 't'==b)
                currentStruct = refStruct.T_rev;
            else
                error('Illegal base');
            end
        else
            error('Exception.');
        end
        
        % Translate/rotation/modify AtomSerNo
        nAtom = length(currentStruct.Model.Atom);
        for j = 1:nAtom
            currentAtom = currentStruct.Model.Atom(j);
            
            % Modify AtomSerNo
            currentAtom.AtomSerNo = firstAtomSerNo_out;
            firstAtomSerNo_out = firstAtomSerNo_out+1;
            
            % Modify resSeq
            %currentAtom.resSeq = tourIndex(i);
            currentAtom.resSeq = i;
            
            % Rotation/translation
            tmpCoor = [currentAtom.X; currentAtom.Y; currentAtom.Z];
            %         if(mod(h,2)==0)
            %             tmpCoor = Ry180*tmpCoor;
            %         end
            tmpCoor = R*tmpCoor + D;
            
            % Append currentAtom to pdbStruct
            currentAtom.X = tmpCoor(1);
            currentAtom.Y = tmpCoor(2);
            currentAtom.Z = tmpCoor(3);
            pdbStruct.Model.Atom(atomIndex) = currentAtom;
            atomIndex = atomIndex+1;
        end
    end
end
pdbStruct.Model.Atom(atomIndex:end) = [];

% Build the terminal for the PDB structure
pdbStruct.Model.Terminal.SerialNo = firstAtomSerNo_out;
firstAtomSerNo_out = firstAtomSerNo_out+1;

pdbStruct.Model.Terminal.resName = currentAtom.resName;
pdbStruct.Model.Terminal.chainID = currentAtom.chainID;
pdbStruct.Model.Terminal.resSeq = currentAtom.resSeq;
pdbStruct.Model.Terminal.iCode = currentAtom.iCode;
end


% Read reference structures
function [pdbStruct_for,pdbStruct_rev] = readReference(filename)
% In the reference structure:
% Running forward: chainID==S
% Running backward: chainID==A

pdbStruct = pdbread(filename);
pdbStruct_for = [];
pdbStruct_for.Model.Atom = [];
pdbStruct_rev = [];
pdbStruct_rev.Model.Atom = [];

% Rotation matrix: first, -90 deg about X-axis; second, -90 deg about
% Y-axis
R = Ry(-90)*Rx(-90);

% Running through all the atoms
nAtom = length(pdbStruct.Model.Atom);
for i = 1:nAtom
    currentAtom = pdbStruct.Model.Atom(i);
    
    % Looking for the forward reference structure
    if('S'==currentAtom.chainID && 2==currentAtom.resSeq)
        currentAtom.AtomSerNo = -1;
        currentAtom.resSeq = -1;
        
        tmpCoor = [currentAtom.X; currentAtom.Y; currentAtom.Z];
        tmpCoor = R*tmpCoor;
        currentAtom.X = tmpCoor(1);
        currentAtom.Y = tmpCoor(2);
        currentAtom.Z = tmpCoor(3);
        
        currentAtom.chainID = 'A';
        pdbStruct_for.Model.Atom = [pdbStruct_for.Model.Atom, currentAtom];
    end
    
    % Looking for the backward reference structure
    if('A'==currentAtom.chainID && 2==currentAtom.resSeq)
        currentAtom.AtomSerNo = -1;
        
        tmpCoor = [currentAtom.X; currentAtom.Y; currentAtom.Z];
        tmpCoor = R*tmpCoor;
        currentAtom.X = tmpCoor(1);
        currentAtom.Y = tmpCoor(2);
        currentAtom.Z = tmpCoor(3);
        
        pdbStruct_rev.Model.Atom = [pdbStruct_rev.Model.Atom, currentAtom];
    end
end

end


% Rotation about X-axis
function R = Rx(deg)
rad = deg/180*pi;
R = [1 0 0; 0 cos(rad) -sin(rad); 0 sin(rad) cos(rad)];
end


% Rotation about Y-axis
function R = Ry(deg)
rad = deg/180*pi;
R = [cos(rad) 0 sin(rad); 0 1 0; -sin(rad) 0 cos(rad)];
end