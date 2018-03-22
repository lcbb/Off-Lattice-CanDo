function [strE] = getStrainEnergyADINAoutTimeStep_v2(outFN,nelem,timestep)

fid = fopen(outFN,'r');
strE = zeros(nelem,1);

% headerLastStep = '  P R I N T   O U T   F O R   T I M E   ( L O A D )  S T E P       450';
% headerStress = '    S T R E S S    C A L C U L A T I O N S';
% headerStrainEnergy = 'STRAIN-ENERGY';
% headerBreak = '1PROGRAM ADINA - VERSION 8.7.3 (build 12.13.2010)   DNA';

headerLastStep = sprintf('S T E P%10d',timestep);
headerStress = 'S T R E S S';
headerStrainEnergy = 'STR-ENERGY';
headerBreak = '1PROGRAM ADINA';

% go to print out for the last time step
while 1
    ss = fgetl(fid);
    if ~isempty(strfind(ss,headerLastStep))
        break;
    end
end

% go to stress part
while 1
    ss = fgetl(fid);
    if ~isempty(strfind(ss,headerStress))
        break;
    end
end

% go to first strain-energy
while 1
    ss = fgetl(fid);
    if ~isempty(strfind(ss,headerStrainEnergy))
        break;
    end
end

% read strain-energy
ielem = 0;
flag = 1;
while flag
    
    % skip blank line
    ss = fgetl(fid);

    % check if data line or not
    ss = fgetl(fid);
    if ~isempty(strfind(ss,headerBreak))
        while 1
            ss = fgetl(fid);
            if ~isempty(strfind(ss,headerStrainEnergy))
                break;
            end
        end
    else
        % read strain energy
        [aaa,bbb] = strread(ss,'%s%*s%*s%*s%*s%*s%*s%*s%*s%s');
        ielem = ielem + 1;
        enum = str2double(char(aaa));
        if isempty(strfind(char(bbb),'E'))
            estr = 0.0;
        else
            estr = str2double(char(bbb));
        end
        
        strE(enum,1) = estr;
        %fprintf(1,'%d\t%f\n',ielem,tmp);
        % skip blank line
        ss = fgetl(fid);
        if ielem == nelem
            flag = 0;
        end    
    end

end


fclose(fid);

