function [totStrE,beamStrE,trussStrE,alignStrE_HJE_ds,alignStrE_HJE_ss,alignStrE_DSDNA,alignStrE_BBE] = getTotalStrE_finalStep(matPATH,outPATH,param)

load(matPATH,'FEModel');

node = FEModel.node(:,1:3);
if(isfield(FEModel,'rigidLinks'))
    elem = [ ...
            FEModel.beamsNormal; ...
            FEModel.beamsNick; ...
            FEModel.rigidLinks; ...
           ];
else
    elem = [ ... 
            FEModel.beamsNormal; ...
            FEModel.beamsNick; ...
            FEModel.beamsBulge; ...
            FEModel.beamsNone; ...
           ];
end
elemBeam = [ ... 
            FEModel.beamsNormal; ...
            FEModel.beamsNick; ...
            FEModel.beamsBulge; ...
            FEModel.beamsNone; ...
           ];   
elemRender = [ ...
               FEModel.beamsNormal; ...
               FEModel.beamsNick; ...
               FEModel.beamsBulge; ...
               FEModel.beamsNone; ...
%                FEModel.shrinkBeamsNormal; ...
%                FEModel.shrinkBeamsNick; ...
             ];
if(isfield(FEModel,'rigidLinks'))
    elemRenderTruss = [ ...
                       FEModel.rigidLinks; ...
                       FEModel.entSprings; ...
                      ];
else
    elemRenderTruss = [FEModel.entSprings];
end
nNode = size(node,1);
nElem = size(elem,1);
nElemBeam = size(elemBeam,1);
nElemRender = size(elemRender,1);
nElemRenderTruss = size(elemRenderTruss,1);
nAlignHJE_ds = numel(find(strcmp(FEModel.typeHJE, 'ds')));
nAlignHJE_ss = numel(find(strcmp(FEModel.typeHJE, 'ss')));
nAlignDSDNA = size(FEModel.alignDSDNA,1);
nAlignBBE = size(FEModel.BBE, 1) + size(FEModel.nickOpen,1);
nTruss = size(FEModel.entSprings,1);


% Retrieve mechanical energy information from the ADINA *.out file
[strE,trussSS,trussEE,misalignHJE_ds,misalignHJE_ss,misalignDSDNA,misalignBBE,reaction,solutionTime_outfile] = ...
    procAdinaOut(outPATH,nElemBeam,nTruss,nAlignHJE_ds,nAlignHJE_ss,nAlignDSDNA,nAlignBBE,param.timeStep);
strE = strE(:,end);
trussSS = trussSS(:,end);
trussEE = trussEE(:,end);
misalignHJE_ds = misalignHJE_ds(:,:,end);
misalignHJE_ss = misalignHJE_ss(:,:,end);
misalignDSDNA = misalignDSDNA(:,:,end);
misalignBBE = misalignBBE(:,:,end);

strE = strE / param.KbT;
beamStrE = sum(strE);
assert(numel(beamStrE)==1);

trussStrE = getTrussStrE(matPATH,param,trussSS,trussEE);
trussStrE = trussStrE / param.KbT;
trussStrE = sum(trussStrE);
assert(numel(trussStrE)==1);

[alignStrE_HJE_ds,alignStrE_HJE_component_ds] = alignE(misalignHJE_ds,param.alignHJE_ds);
alignStrE_HJE_ds = alignStrE_HJE_ds/param.KbT;
alignStrE_HJE_component_ds = alignStrE_HJE_component_ds/param.KbT;
alignStrE_HJE_ss = alignE(misalignHJE_ss,param.alignHJE_ss)/param.KbT;
alignStrE_DSDNA = alignE(misalignDSDNA,param.alignDSDNA)/param.KbT;
alignStrE_BBE = sum(misalignBBE(2,:)*param.alignBBE)/param.KbT;

totStrE = beamStrE + trussStrE + alignStrE_HJE_ds + alignStrE_HJE_ss + alignStrE_DSDNA + alignStrE_BBE;
% totStrE = beamStrE + trussStrE + sum(alignStrE_HJE_component(4:6));
assert(numel(totStrE)==1);

end


% Mechanical energy for alignment elements
function [E,E_component] = alignE(misalign,K)

dof = zeros(6,size(misalign,2));
dof(1:3,:) = misalign(2:4,:);
dof(4:6,:) = bsxfun(@times,misalign(5,:),misalign(6:8,:));
E_component(1) = sum(1/2 * (K.KTC1*dof(1,:).^2));
E_component(2) = sum(1/2 * (K.KTC2*dof(2,:).^2));
E_component(3) = sum(1/2 * (K.KTC3*dof(3,:).^2));
E_component(4) = sum(1/2 * (K.KRC1*dof(4,:).^2));
E_component(5) = sum(1/2 * (K.KRC2*dof(5,:).^2));
E_component(6) = sum(1/2 * (K.KRC3*dof(6,:).^2));
% 1/2 * (K.KTC1*dof(1,:).^2 + K.KTC2*dof(2,:).^2 + K.KTC3*dof(3,:).^2 + ...
%                      K.KRC1*dof(4,:).^2 + K.KRC2*dof(5,:).^2 + K.KRC3*dof(6,:).^2);
E = sum(E_component);

end