function trussStrE = getTrussStrE(matPATH,param,trussSS,trussEE)

% Load the FE model
load(matPATH, 'FEModel','connSSE');

nTruss = size(FEModel.entSprings,1);
assert(nTruss==size(connSSE,1) && nTruss==size(trussSS,1) && nTruss==size(trussEE,1));
nTimeStep = size(trussSS,2);
assert(nTimeStep==size(trussEE,2));


% Nonlinear stress-strain relationship (See runFE_offlattice.m, Ln 392-398)
ss1 = -10:9:-1;                         % Stress
ss2 = [1:10 15:5:50 60:10:500 1e20];    % Stress
ee1 = zeros(nTruss,numel(ss1));         % Strain
ee2 = zeros(nTruss,numel(ss2));         % Strain
R0 = zeros(nTruss,1);                   % Initial length
N0 = zeros(nTruss,1);                   % Number of bases
for i = 1:nTruss
    nid1 = FEModel.entSprings(i,2);
    nid2 = FEModel.entSprings(i,3);
    nBase = connSSE(i,5);
    distNode = norm(FEModel.node(nid1,1:3) - FEModel.node(nid2,1:3));
    if(nBase >= 5)
        ee1(i,:) = (5*nBase)/distNode*(coth(ss1*15/param.KbT)-param.KbT./ss1/15).*(1+ss1/800) - 1.0;
        ee2(i,:) = (5*nBase)/distNode*(coth(ss2*15/param.KbT)-param.KbT./ss2/15).*(1+ss2/800) - 1.0;
    else
        ee1(i,:) = (3.4*nBase+3.4)/distNode.*(1+ss1/800) - 1.0;
        ee2(i,:) = (3.4*nBase+3.4)/distNode.*(1+ss2/800) - 1.0;
    end
    R0(i) = distNode;
    N0(i) = nBase;
end
ee = [-ones(nTruss,1), ee2];
ss = [0.0, ss2];


% Calculate strain energy by integration
trussStrE = zeros(nTruss,nTimeStep);
for i = 1:nTimeStep
    trussStrE(:,i) = getTrussStrE_singleStep(trussSS(:,i),trussEE(:,i),ss,ee,R0,N0);
end

end



% Calculate strain energy for a single time step
function trussStrE = getTrussStrE_singleStep(trussSS,trussEE,ss,ee,R0,N0)

assert(iscolumn(trussSS) && iscolumn(trussEE));
nTruss = size(trussSS,1);

trussStrE = zeros(nTruss,1);
for i = 1:nTruss
    if(N0(i) >= 5)
        tmp1 = sum(ss<=trussSS(i));
        tmp2 = sum(ee(i,:)<=trussEE(i));
        assert(tmp1==tmp2 && tmp1>=1 && tmp1<size(ss,2));
        currRange = tmp1;
        
        for j = 1:(currRange-1)
            ss_left = ss(j);
            ss_right = ss(j+1);
            ee_left = ee(i,j);
            ee_right = ee(i,j+1);
            trussStrE(i) = trussStrE(i) + R0(i) * 1/2*(ss_left+ss_right)*(ee_right-ee_left);
        end
        trussStrE(i) = trussStrE(i) + R0(i) * 1/2*(ss(currRange)+trussSS(i))*(trussEE(i)-ee(currRange));
    else
        trussStrE(i) = R0(i) * 1/2*trussSS(i)*trussEE(i);
    end
    assert(trussStrE(i) >= -1e-6);
end

end