function triad = getTriad(matPATH,porPATH,solutionTime)

rot = [0 1 0;0 0 1;1 0 0];  % Convert CanDo convention to 3DNA convention
initTriad = getInitialTriad(matPATH);
nNode = size(initTriad,1);

if(isempty(solutionTime))
    % solutionTime = 0
    triad = cell(size(initTriad));
    for i = 1:nNode
        triad{i} = initTriad{i} * rot;
    end
else
    % solutionTime > 0
    currTriad = getDeformedTriad(matPATH,porPATH,solutionTime);
    nSolutionTime = size(currTriad,2);
    assert(nNode == size(currTriad,1));
    
    triad = cell(size(currTriad));
    for loop = 1:nSolutionTime
        for i = 1:nNode
            triad{i,loop} = currTriad{i,loop} * initTriad{i} * rot;
        end
    end
end

end


function triad = getInitialTriad(matPATH)

load(matPATH,'FEModel');
nNode = size(FEModel.node,1);

triad = cell(nNode,1);
for i = 1:nNode
    e1 = FEModel.nodalTriad(i,:,1);
    e1 = e1/norm(e1);
    e2 = FEModel.nodalTriad(i,:,2);
    e2 = e2/norm(e2);
    e3 = cross(e1,e2);
    triad{i} = [e3', e1', e2'];
end

end


function triad = getDeformedTriad(matPATH,porPATH,solutionTime)

nSolutionTime = numel(solutionTime);
% Read in the finite element model
load(matPATH,'FEModel');
elemBeam = [ ... 
             FEModel.beamsNormal; ...
             FEModel.beamsNick; ...
           ];
nElem = size(elemBeam,1);
nNode = size(FEModel.node,1);

% Get quaternions
[initriad,curtriad] = getTriadQuaternion(FEModel,porPATH,solutionTime);

% Convert quaternions to triads
triad = cell(nNode,nSolutionTime);
for loop = 1:nSolutionTime
    % Get triads for all the elements
    elemTriad_init = quat2dcm(cell2mat(initriad));
    elemTriad1_curr = quat2dcm(cell2mat(curtriad(:,1,loop)));
    elemTriad2_curr = quat2dcm(cell2mat(curtriad(:,2,loop)));
    for i = 1:nElem
        iNode1 = elemBeam(i,2);
        iNode2 = elemBeam(i,3);
        elemTriad1 = elemTriad1_curr(:,:,i) * elemTriad_init(:,:,i)';
        elemTriad2 = elemTriad2_curr(:,:,i) * elemTriad_init(:,:,i)';
        
        if(isempty(triad{iNode1,loop}))
            triad{iNode1,loop} = elemTriad1;
        else
%             assert(norm(triad{iNode1,loop}-elemTriad1, 'fro') < 1e-6);
        end
        
        if(isempty(triad{iNode2,loop}))
            triad{iNode2,loop} = elemTriad2;
        else
%             assert(norm(triad{iNode2,loop}-elemTriad2, 'fro') < 1e-6);
        end
    end
    
    % Verification
    for i = 1:nNode
        assert(~isempty(triad{i,loop}));
    end
end

end


function [initriad,curtriad] = getTriadQuaternion(FEModel,porPATH,solutionTime)

fid = fopen(porPATH);
nSolutionTime = numel(solutionTime);
elemBeam = [ ... 
             FEModel.beamsNormal; ...
             FEModel.beamsNick; ...
           ];
nElem = size(elemBeam,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in INITRIAD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Go to the first keyword 'ELEMENT4'
while(1)
    ss = fgetl(fid);
    if(~isempty(strfind(ss,'ELEMENT4')))
        break;
    end
end

% Read triads
initriad = cell(nElem,1);
for i = 1:nElem
    while(1)
        ss = fgetl(fid);
        if(~isempty(strfind(ss,'INITRIAD')))
            break;
        end
    end
    ss = fgetl(fid);
    initriad{i} = [str2double(ss(1:20)) str2double(ss(21:40)) str2double(ss(41:60)) str2double(ss(61:80))];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in CURTRIAD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
curtriad = cell(nElem,2,nSolutionTime);
for loop = 1:nSolutionTime
    currentSolutionTime = solutionTime(loop);
    
    % Go to the first keyword 'NEWSTEP4'
    while(1)
        ss = fgetl(fid);
        if(~isempty(strfind(ss,'NEWSTEP4')))
            for i = 1:6
                ss = fgetl(fid);
            end
            ctime = str2double(ss(41:60));
            if(abs(ctime-currentSolutionTime)/currentSolutionTime<1.0e-4);
                break;
            end
        end
    end
    
    % Read triads
    for i = 1:nElem
        while(1)
            ss = fgetl(fid);
            if(~isempty(strfind(ss,'CURTRIAD')))
                break;
            end
        end
        fgetl(fid);
        ss = fgetl(fid);
        curtriad{i,1,loop} = [str2double(ss(1:20)) str2double(ss(21:40)) str2double(ss(41:60)) str2double(ss(61:80))];
        ss = fgetl(fid);
        curtriad{i,2,loop} = [str2double(ss(1:20)) str2double(ss(21:40)) str2double(ss(41:60)) str2double(ss(61:80))];
    end
end

end