% Get the coordinates of each node in the solution shape
%     - dNode: the actual Cartesian coordinates in the solution shape
%     - pNode: the Cartesian coordinates of the solution shape projected
%              onto the principal axes
%
function [dNode,pNode,pAxes] = getSolutionShape(matPATH,displ)

load(matPATH,'FEModel');
dNode = FEModel.node(:,1:3) + displ(:,1:3);
[pAxes,pNode] = pca(dNode);
if(dot(cross(pAxes(:,1),pAxes(:,2)),pAxes(:,3))<0);
    pAxes(:,3) = -pAxes(:,3);
    pNode(:,3) = -pNode(:,3);
end

end