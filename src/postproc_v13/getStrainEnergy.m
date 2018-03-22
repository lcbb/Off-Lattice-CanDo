function strE = getStrainEnergy(outPATH,matPATH,timestep,param)

load(matPATH,'FEModel');
nelem = size(FEModel.beamsNormal,1) + size(FEModel.beamsNick,1);
strE = getStrainEnergyADINAoutTimeStep_v2(outPATH,nelem,timestep);
strE = strE / param.KbT;

end