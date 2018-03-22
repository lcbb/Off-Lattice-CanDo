function [eVal] = getEigenvaluesADINA(feFNout,startMode,endMode)

nMode = endMode-startMode+1;

eVal = zeros(nMode,1);
fid = fopen(feFNout,'r');
while 1
    ss = fgetl(fid);
    if strfind(ss,'FREQUENCY (RAD/SEC)      FREQUENCY (CYCLES/SEC)        PERIOD (SECONDS)')
        break;
    end
end
for itmp = 1:startMode
    ss = fgetl(fid);
end
for itmp = 1:nMode
    ss = fgetl(fid);
    tmp = sscanf(ss,'%f');
    eVal(itmp,1) = tmp(2)^2;  % tmp(2) = frequency (rad/sec)
end
fclose(fid);