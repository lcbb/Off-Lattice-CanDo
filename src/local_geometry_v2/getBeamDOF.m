function beamDOF = getBeamDOF(FEModelPATH,dNode,triad)

% Read in the finite element model
load(FEModelPATH,'FEModel');
elemBeam = [ ... 
             FEModel.beamsNormal; ...
             FEModel.beamsNick; ...
           ];
nElem = size(elemBeam,1);
nNode = size(FEModel.node,1);

assert(nNode==size(dNode,1) && 3==size(dNode,2));
assert(nNode==numel(triad));

% Calculate the 6 DOFs for a beam element using the 3DNA convention
beamDOF.shift = zeros(nElem,1);
beamDOF.slide = zeros(nElem,1);
beamDOF.rise  = zeros(nElem,1);
beamDOF.tilt  = zeros(nElem,1);
beamDOF.roll  = zeros(nElem,1);
beamDOF.twist = zeros(nElem,1);

for i = 1:nElem
    iNode1 = elemBeam(i,2);
    iNode2 = elemBeam(i,3);
    
    if(dot(triad{iNode1}(:,3), dNode(iNode2,:)'-dNode(iNode1,:)') > 0)
        [beamDOF.shift(i),beamDOF.slide(i),beamDOF.rise(i), ...
         beamDOF.tilt(i),beamDOF.roll(i),beamDOF.twist(i)] = ...
         baseStepParam_v2(triad{iNode1},dNode(iNode1,:)',triad{iNode2},dNode(iNode2,:)','basepairstep');
    end
end

end