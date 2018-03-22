function [] = generateBILD(dnaTop, strand, conn_dsDNA_ssDNA, HJE, filename, option)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R_base = 1.5;
R_backbone = 0.5;
R_basepair = 0.25;
R_triad = 0.5;

R_cone = 1.8;
L_cone = 3.6;

colorList = {'tan', 'salmon', 'orange', 'gold', 'sea green', 'steel blue', 'medium purple', 'rosy brown'};
defaultColor = 'light gray';
alignColor = 'firebrick';
connColor = 'firebrick';

% Colors for four bases
% Reference: http://www.umass.edu/molvis/tutorials/dna/atgc.htm
AzureRGB = [0, 127, 255]/255;       % color for A
TweetyBirdRGB = [253, 245, 1]/255;  % color for T
GreenRGB = [0, 255, 0]/255;         % color for G
CarmineRGB = [150, 0, 24]/255;      % color for C

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write the file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(filename, 'w');

% Draw the structure
for i = 1:numel(strand)
    lenTour = numel(strand(i).tour);

    for j = 1:lenTour
        % Draw the base
        id = strand(i).tour(j);
        if(strcmp(option, 'strand'))
            fprintf(fid, '.color %s\n', colorList{mod(i-1,numel(colorList))+1});
        elseif(strcmp(option, 'HJE'))
            k = findHJE(id, HJE);
            if(numel(k) == 1)
                fprintf(fid, '.color %s\n', colorList{mod(k-1,numel(colorList))+1});
            elseif(numel(k) == 2)
                fprintf(fid, '.color %s\n', alignColor);
            else
                fprintf(fid, '.color %s\n', defaultColor);
            end
        elseif(strcmp(option, 'sequence'))
            currSeq = dnaTop(id).seq;
            if(strcmp(currSeq, 'A'))
                fprintf(fid, '.color %f %f %f\n', AzureRGB(1), AzureRGB(2), AzureRGB(3));
            elseif(strcmp(currSeq, 'T'))
                fprintf(fid, '.color %f %f %f\n', TweetyBirdRGB(1), TweetyBirdRGB(2), TweetyBirdRGB(3));
            elseif(strcmp(currSeq, 'G'))
                fprintf(fid, '.color %f %f %f\n', GreenRGB(1), GreenRGB(2), GreenRGB(3));
            elseif(strcmp(currSeq, 'C'))
                fprintf(fid, '.color %f %f %f\n', CarmineRGB(1), CarmineRGB(2), CarmineRGB(3));
            else
                fprintf(fid, '.color %s\n', defaultColor);
            end
        elseif(strcmp(option, 'conn_dsDNA_ssDNA'))
            fprintf(fid, '.color %s\n', defaultColor);
        else
            error('Illegal option');
        end
        fprintf(fid, '.sphere\t%.3f\t%.3f\t%.3f\t%.3f\n', ...
                dnaTop(id).xyz(1)*10, dnaTop(id).xyz(2)*10, dnaTop(id).xyz(3)*10, R_base);
            
        % Draw the backbones
        idDown = dnaTop(id).down;
        if(idDown >= 0)
            if(strcmp(option, 'conn_dsDNA_ssDNA') && ~isempty(find(conn_dsDNA_ssDNA(:,4)==id & conn_dsDNA_ssDNA(:,8)==idDown,1)))
                fprintf(fid, '.color %s\n', connColor);
            else
                fprintf(fid, '.color %s\n', defaultColor);
            end
            fprintf(fid, '.cylinder\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n', ...
                    dnaTop(id).xyz(1)*10, dnaTop(id).xyz(2)*10, dnaTop(id).xyz(3)*10, ...
                    dnaTop(idDown).xyz(1)*10, dnaTop(idDown).xyz(2)*10, dnaTop(idDown).xyz(3)*10, ...
                    R_backbone);
        end
        
        % Draw the Watson-Crick connections
        idAcross = dnaTop(id).across;
        if(idAcross >= 0)
            fprintf(fid, '.color %s\n', defaultColor);
            fprintf(fid, '.cylinder\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n', ...
                    dnaTop(id).xyz(1)*10, dnaTop(id).xyz(2)*10, dnaTop(id).xyz(3)*10, ...
                    dnaTop(idAcross).xyz(1)*10, dnaTop(idAcross).xyz(2)*10, dnaTop(idAcross).xyz(3)*10, ...
                    R_basepair);
        end
    end
    
    % Draw the strand directionality
    if(lenTour > 1)
        id = strand(i).tour(lenTour);
        if(strand(i).isCircular)
            idDown = strand(i).tour(1);
            V_cone = dnaTop(idDown).xyz - dnaTop(id).xyz;
        else
            idUp = strand(i).tour(lenTour-1);
            V_cone = dnaTop(id).xyz - dnaTop(idUp).xyz;
        end
        V_cone = V_cone/norm(V_cone);
    else
        V_cone = [1,0,0]';
    end
    vertex_cone = dnaTop(id).xyz*10 + V_cone*L_cone;
    
    if(strcmp(option, 'strand'))
        fprintf(fid, '.color %s\n', colorList{mod(i-1,numel(colorList))+1});
    elseif(strcmp(option, 'HJE'))
        k = findHJE(id, HJE);
        if(numel(k) == 1)
            fprintf(fid, '.color %s\n', colorList{mod(k-1,numel(colorList))+1});
        elseif(numel(k) == 2)
            fprintf(fid, '.color %s\n', alignColor);
        else
            fprintf(fid, '.color %s\n', defaultColor);
        end
    elseif(strcmp(option, 'sequence'))
        currSeq = dnaTop(id).seq;
        if(strcmp(currSeq, 'A'))
            fprintf(fid, '.color %f %f %f\n', AzureRGB(1), AzureRGB(2), AzureRGB(3));
        elseif(strcmp(currSeq, 'T'))
            fprintf(fid, '.color %f %f %f\n', TweetyBirdRGB(1), TweetyBirdRGB(2), TweetyBirdRGB(3));
        elseif(strcmp(currSeq, 'G'))
            fprintf(fid, '.color %f %f %f\n', GreenRGB(1), GreenRGB(2), GreenRGB(3));
        elseif(strcmp(currSeq, 'C'))
            fprintf(fid, '.color %f %f %f\n', CarmineRGB(1), CarmineRGB(2), CarmineRGB(3));
        else
            fprintf(fid, '.color %s\n', defaultColor);
        end
    elseif(strcmp(option, 'conn_dsDNA_ssDNA'))
        fprintf(fid, '.color %s\n', defaultColor);
    else
        error('Illegal option');
    end
    fprintf(fid, '.cone\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n', ...
            dnaTop(id).xyz(1)*10, dnaTop(id).xyz(2)*10, dnaTop(id).xyz(3)*10, ...
            vertex_cone(1), vertex_cone(2), vertex_cone(3), ...
            R_cone);
end

% Draw the reference frame for each element
for i = 1:numel(HJE)
    c = HJE(i).c;
    e1_end = HJE(i).c + HJE(i).e1;
    e2_end = HJE(i).c + HJE(i).e2;
    e3_end = HJE(i).c + HJE(i).e3;
    
    fprintf(fid, '.color red\n');
    fprintf(fid, '.arrow\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n', ...
            c(1)*10, c(2)*10, c(3)*10, e1_end(1)*10, e1_end(2)*10, e1_end(3)*10, R_triad, R_triad*2);

    fprintf(fid, '.color green\n');
    fprintf(fid, '.arrow\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n', ...
            c(1)*10, c(2)*10, c(3)*10, e2_end(1)*10, e2_end(2)*10, e2_end(3)*10, R_triad, R_triad*2);

    fprintf(fid, '.color blue\n');
    fprintf(fid, '.arrow\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n', ...
            c(1)*10, c(2)*10, c(3)*10, e3_end(1)*10, e3_end(2)*10, e3_end(3)*10, R_triad, R_triad*2);
end

fclose(fid);

end


function k = findHJE(id, HJE)

k = [];
for i = 1:numel(HJE)
    tmp = ([HJE(i).idSeq_main; HJE(i).idSeq_comp] == id);
    k = cat(1, k, i*ones(sum(tmp),1));
end
assert(numel(k) <= 2);

end