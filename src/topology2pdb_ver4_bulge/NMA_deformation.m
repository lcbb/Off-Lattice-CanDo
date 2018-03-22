% function [R, D] = NMA_deformation(CanDoPath, filename, T, nFrame, nPeriod)
function [R, D] = NMA_deformation(CanDoPath, filename)
% R - Rotation matrix
% D - Translation
%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Constants and parameters
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kBT = 1.3806504e-1 * T;
% 
% startMode = 1;      % main_CanDoWebServerVer2_proc1.m, Ln 369
% endMode = 206;      % main_CanDoWebServerVer2_proc1.m, Ln 370
% rigid_mode_idx = [1 2 3 4 5 6];     % renderNMABildFormat.m, Ln 44
% 
% load(fullfile(CanDoPath, [filename, '_FEModel.mat']));
% % generateAdinaIn4NMA.m, Ln 14-16
% node = FEModel.node(:,1:3);
% nNode = size(node,1);
% endMode = min(endMode,nNode*6);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Eigenvalues and eigenvectors
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eVal = getEigenvalue(CanDoPath, filename, startMode, endMode);
% eVec = getEigenvector(CanDoPath, filename, nNode, startMode, endMode);
% 
% nMode = endMode-startMode+1-numel(rigid_mode_idx);
% eVal(rigid_mode_idx,:) = [];
% eVec(:,:,rigid_mode_idx) = [];
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Generate the frames
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% phi = zeros(nNode,6,nFrame);
% tmax = nPeriod*(2*pi/sqrt(eVal(1,1)));
% dt = tmax/nFrame;
% t = dt:dt:tmax;
% 
% p0 = 2*pi*rand(nMode,1);
% for i = 1:nFrame
%     for j = 1:nMode
%         phi(:,:,i) = phi(:,:,i) + sqrt(2*kBT/eVal(j))*eVec(:,:,j)*cos(sqrt(eVal(j,1))*t(i)+p0(j));
%     end
% end

load(fullfile(CanDoPath,strcat(filename,'_phi.mat')), 'phi');
nNode = size(phi,1);
assert(6 == size(phi,2));
nFrame = size(phi,3);

R = zeros(3,3,nNode,nFrame);    % Rotation matrix
D = zeros(3,nNode,nFrame);    % Translation vector
for i = 1:nFrame
    for j = 1:nNode
        tmp = zeros(1,4);
        tmp(1:3) = phi(j,4:6,i);
        tmp(4) = norm(phi(j,4:6,i));
        tmp(1:3) = tmp(1:3)/tmp(4);
        R(:,:,j,i) = vrrotvec2mat(tmp);
        D(:,j,i) = phi(j,1:3,i)';
    end
end

end


% % Modified from getEigenvaluesADINA.m
% function eVal = getEigenvalue(CanDoPath,filename,startMode,endMode)
% 
% feFNout = fullfile(CanDoPath, [filename, '_NMA.out']);
% 
% nMode = endMode-startMode+1;
% 
% eVal = zeros(nMode,1);
% fid = fopen(feFNout,'r');
% while 1
%     ss = fgetl(fid);
%     if strfind(ss,'FREQUENCY (RAD/SEC)      FREQUENCY (CYCLES/SEC)        PERIOD (SECONDS)')
%         break;
%     end
% end
% for itmp = 1:startMode
%     ss = fgetl(fid);
% end
% for itmp = 1:nMode
%     ss = fgetl(fid);
%     tmp = sscanf(ss,'%f');
%     eVal(itmp,1) = tmp(2)^2;  % tmp(2) = frequency (rad/sec)
% end
% fclose(fid);
% 
% end
% 
% 
% % Modified from getEigenvectorsADINA.m
% function [eVec] = getEigenvector(CanDoPath,filename,nNode,startMode,endMode)
% 
% listFN = fullfile(CanDoPath, [filename, '_NMA_evec.list']);
% 
% % read eigenvectors
% fid = fopen(listFN,'r');
% eVec = zeros(nNode,6,endMode-startMode+1);
% for i=1:endMode-startMode+1
%     while 1
%         ss = fgetl(fid);
%         if strfind(ss,'Mode_number')
%             break;
%         end
%     end
%     ss = fgetl(fid);
%     for j=1:nNode
%         ss = fgetl(fid);
%         eVec(j,:,i) = cell2mat(textscan(ss,'%*s%*s%f%f%f%f%f%f'));
%     end
% end
% fclose(fid);
% 
% end