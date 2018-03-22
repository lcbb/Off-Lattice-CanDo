function [strE,trussSS,trussEE,misalignHJE_ds,misalignHJE_ss,misalignDSDNA,misalignBBE,reaction,solutionTime] = procAdinaOut(outFN,nBeam,nTruss,nAlignHJE_ds,nAlignHJE_ss,nAlignDSDNA,nAlignBBE,timeStep)

% Headers in the *.out file
headerStep         = '  P R I N T   O U T   F O R   T I M E   ( L O A D )  S T E P';
headerStress       = '    S T R E S S    C A L C U L A T I O N S';
headerStrainEnergy = 'STR-ENERGY';
headerBreak        = '1PROGRAM ADINA';
headerStressTruss  = ' S T R E S S   C A L C U L A T I O N S   F O R   E L E M E N T   G R O U P%7d  ( TRUSS )';
headerAlignHJE_ds  = sprintf(' R E S U L T S   F O R   E L E M E N T   G R O U P%10d   ( ALIGNMENT )', 4+nTruss+1);
headerAlignHJE_ss  = sprintf(' R E S U L T S   F O R   E L E M E N T   G R O U P%10d   ( ALIGNMENT )', 4+nTruss+2);
headerAlignDSDNA   = sprintf(' R E S U L T S   F O R   E L E M E N T   G R O U P%10d   ( ALIGNMENT )', 4+nTruss+3);
headerAlignBBE     = sprintf(' R E S U L T S   F O R   E L E M E N T   G R O U P%10d   ( ALIGNMENT )', 4+nTruss+4);
headerReaction     = '   P R I N T O U T   O F   R E A C T I O N S';
headerRestart      = '  R E S T A R T    D A T A    S A V I N G    I N F O R M A T I O N';

% Initialization
nTimeStep = numel(timeStep);
strE = zeros(nBeam,nTimeStep);
trussSS = zeros(nTruss,nTimeStep);
trussEE = zeros(nTruss,nTimeStep);
misalignHJE_ds = zeros(8,nAlignHJE_ds,nTimeStep);
misalignHJE_ss = zeros(8,nAlignHJE_ss,nTimeStep);
misalignDSDNA = zeros(8,nAlignDSDNA,nTimeStep);
misalignBBE = zeros(2,nAlignBBE,nTimeStep);
reaction = [];
solutionTime = zeros(nTimeStep,1);

fid = fopen(outFN,'r');
for i = 1:nTimeStep
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Go to the current time step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    curr_headerStep = sprintf('%s%10d', headerStep,timeStep(i));
    while(1)
        ss = fgetl(fid);
        if(~ischar(ss))
            error('Cannot find string\n%s\n', curr_headerStep);
        end
        if(~isempty(strfind(ss,curr_headerStep)))
            break;
        end
    end
%     fprintf('Reached time step %d\n', i);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Go to the beam stress calculation section in the current time step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while(1)
        ss = fgetl(fid);
        if(~ischar(ss))
            error('Cannot find string\n%s\n', headerStress);
        end
        if(~isempty(strfind(ss,headerStress)))
            break;
        end
    end
%     fprintf('Reached stress calculation section\n');
    
    % Go to the first beam strain energy
    while(1)
        ss = fgetl(fid);
        if(~ischar(ss))
            error('Cannot find string\n%s\n', headerStrainEnergy);
        end
        if(~isempty(strfind(ss,headerStrainEnergy)))
            break;
        end
    end
%     fprintf('Reached the first beam strain energy\n');
    
    % Read strain energy
    nRead = 0;
    while(nRead<nBeam)
        % Skip one line
        fgetl(fid);
        
        % Check if data line or not
        ss = fgetl(fid);
        if(~isempty(strfind(ss,headerBreak)))
            % Find a break, and go to the next strain energy header
            while(1)
                ss = fgetl(fid);
                if(~ischar(ss))
                    error('Cannot find string\n%s\n', headerStrainEnergy);
                end
                if ~isempty(strfind(ss,headerStrainEnergy))
                    break;
                end
            end
            
        else
            % Read strain energy
            nRead = nRead+1;
            raw = textscan(ss,'%s%*s%*s%*s%*s%*s%*s%*s%*s%s');
            enum = str2double(raw{1});
            if isempty(strfind(raw{2},'E'))
                estr = 0.0;
            else
                estr = str2double(raw{2});
            end
            strE(enum,i) = estr;
            
            % Skip one line
            fgetl(fid);
        end
    end
%     fprintf('Read all strain energy\n');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Go to the truss stress calculation section in the current time step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    for j = 1:nTruss
        % Header for the current truss element
        headerStressTruss_curr = sprintf(headerStressTruss, 4+j);
        
        % Go to the current truss element
        while(1)
            ss = fgetl(fid);
            if(~ischar(ss))
                error('Cannot find string\n%s\n', headerStressTruss_curr);
            end
            if(~isempty(strfind(ss,headerStressTruss_curr)))
                break;
            end
        end
        
        % Skip 9 lines
        for k = 1:9
            fgetl(fid);
        end
        
        % Read truss stress/strain
        ss = fgetl(fid);
        raw = textscan(ss,'%s%s%s%s');
        trussSS(j,i) = str2double(raw{3});
        trussEE(j,i) = str2double(raw{4});
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Go to the misalignment section in the current time step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Alignment elements for HJE
    if(nAlignHJE_ds>0)
        % Go to the first misalignment
        while(1)
            ss = fgetl(fid);
            if(~ischar(ss))
                error('Cannot find string\n%s\n', headerAlignHJE_ds);
            end
            if(~isempty(strfind(ss,headerAlignHJE_ds)))
                break;
            end
        end
        
        % Skip 11 lines
        for j = 1:11
            fgetl(fid);
        end
        
        % Read misalignments
        for j = 1:nAlignHJE_ds
            fgetl(fid);
            ss = fgetl(fid);
            raw = textscan(ss,'%s%s%s%s%s%s%s%s');
            for k = 1:8
                misalignHJE_ds(k,j,i) = str2double(raw{k});
            end
        end
    end
    
    % Alignment elements for single-stranded crossovers
    if(nAlignHJE_ss>0)
        % Go to the first misalignment
        while(1)
            ss = fgetl(fid);
            if(~ischar(ss))
                error('Cannot find string\n%s\n', headerAlignHJE_ss);
            end
            if(~isempty(strfind(ss,headerAlignHJE_ss)))
                break;
            end
        end
        
        % Skip 11 lines
        for j = 1:11
            fgetl(fid);
        end
        
        % Read misalignments
        for j = 1:nAlignHJE_ss
            fgetl(fid);
            ss = fgetl(fid);
            raw = textscan(ss,'%s%s%s%s%s%s%s%s');
            for k = 1:8
                misalignHJE_ss(k,j,i) = str2double(raw{k});
            end
        end
    end
    
    % Alignment elements for DNA duplex
    if(nAlignDSDNA>0)
        % Go to the first misalignment
        while(1)
            ss = fgetl(fid);
            if(~ischar(ss))
                error('Cannot find string\n%s\n', headerAlignDSDNA);
            end
            if(~isempty(strfind(ss,headerAlignDSDNA)))
                break;
            end
        end
        
        % Skip 11 lines
        for j = 1:11
            fgetl(fid);
        end
        
        % Read misalignments
        for j = 1:nAlignDSDNA
            fgetl(fid);
            ss = fgetl(fid);
            raw = textscan(ss,'%s%s%s%s%s%s%s%s');
            for k = 1:8
                if(~isempty(raw{k}))
                    misalignDSDNA(k,j,i) = str2double(raw{k});
                else
                    misalignDSDNA(k,j,i) = NaN;
                end
            end
        end
    end
    
    % Alignment elements for phosphodiester backbones
%     if(nAlignBBE>0)
%         % Go to the first misalignment
%         while(1)
%             ss = fgetl(fid);
%             if(~ischar(ss))
%                 error('Cannot find string\n%s\n', headerAlignBBE);
%             end
%             if(~isempty(strfind(ss,headerAlignBBE)))
%                 break;
%             end
%         end
%         
%         % Skip 11 lines
%         for j = 1:11
%             fgetl(fid);
%         end
%         
%         % Read misalignments
%         for j = 1:nAlignBBE
%             fgetl(fid);
%             ss = fgetl(fid);
%             raw = textscan(ss,'%s%s');
%             for k = 1:2
%                 misalignBBE(k,j,i) = str2double(raw{k});
%             end
%         end
%     end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Go to the reaction section
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while(1)
        ss = fgetl(fid);
        if(~ischar(ss))
            error('Cannot find string\n%s\n', headerReaction);
        end
        if(~isempty(strfind(ss,headerReaction)))
            break;
        end
    end
%     fprintf('Reached reaction section\n');
    
    % Skip 6 lines
    for j = 1:6
        fgetl(fid);
    end
    
    % Read reaction
    ss = fgetl(fid);
    nnum = 0;
    while(~isempty(ss))
        nnum = nnum+1;
        raw = textscan(ss,'%s%s%s%s%s%s%s%s%s');
        assert(strcmp(raw{2},'GLOBAL') && strcmp(raw{3},'GLOBAL'));
        reaction(nnum,:,i) = [str2double(raw{1}), ...
                            str2double(raw{4}), str2double(raw{5}), str2double(raw{6}), ...
                            str2double(raw{7}), str2double(raw{8}), str2double(raw{9})];
        ss = fgetl(fid);
    end
%     fprintf('Read all reaction\n');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Go to the restart section
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while(1)
        ss = fgetl(fid);
        if(~ischar(ss))
            error('Cannot find string\n%s\n', headerRestart);
        end
        if(~isempty(strfind(ss,headerRestart)))
            break;
        end
    end
%     fprintf('Reached restart section\n');
    
    % Skip one line
    fgetl(fid);
    
    % Get current simulation time
    ss = fgetl(fid);
    raw = regexp(ss,'AT TIME EQUALS','split');
    solutionTime(i) = str2double(raw{2});
%     fprintf('Read simulation time\n');
end
fclose(fid);

end