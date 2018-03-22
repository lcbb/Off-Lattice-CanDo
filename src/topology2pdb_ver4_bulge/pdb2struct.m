function pdbStruct = pdb2struct(pdbName)
% The standard PDB file format is at
% http://deposit.rcsb.org/adit/docs/pdb_atom_format.html

% Initialize the output structure
bufferSize = 1e6;
pdbStruct = [];
pdbStruct.chainSerNo = zeros(bufferSize,1);     % Chain serial number
pdbStruct.chainID = cell(bufferSize,1);         % Chain identifier
pdbStruct.resSeq = zeros(bufferSize,1);         % Residue sequence number
pdbStruct.resName = cell(bufferSize,1);         % Residue name
pdbStruct.AtomSerNo = zeros(bufferSize,1);      % Atom serial number
pdbStruct.AtomName = cell(bufferSize,1);        % Atom name
pdbStruct.XYZ = zeros(bufferSize,3);            % Coordinates

% Read the file line-by-line
fid = fopen(pdbName);

i = 0;
currChainSerNo = 0;
prevChainID = '';
tline = fgetl(fid);
currModel = 1;
while(ischar(tline))
    if(length(tline)>=5 && strcmp(tline(1:5),'MODEL'))
        currModel = str2num(tline(6:end));
    end
    if(length(tline)>=5 && strcmp(tline(1:5),'ATOM '))
        i = i+1;
        
        % 1.1 Read chain serial number
        currChainID = strtrim(tline(21:22));
        if(~strcmp(prevChainID,currChainID))
            currChainSerNo = currChainSerNo+1;
            prevChainID = currChainID;
        end
        %pdbStruct.chainSerNo(i) = currChainSerNo;
        pdbStruct.chainSerNo(i) = currModel;
        % 1.2 Read chain identifier
        pdbStruct.chainID{i} = currChainID;
        
        % 2.1 Read residue sequence number
        pdbStruct.resSeq(i) = str2double(tline(23:26));
        % 2.2 Read residue name, which could be
        % [A, G, C, T], or [DA, DG, DC, DT], or [ADE, GUA, CYT, THY]
        pdbStruct.resName{i} = strtrim(tline(18:20));
        
        % 3.1 Read atom serial number,
        pdbStruct.AtomSerNo(i) = str2double(tline(6:11));
        % 3.2 Read atom name
        pdbStruct.AtomName{i} = strtrim(tline(13:16));
        
        % 4. Read orthogonal coordinates for X/Y/Z in Angstroms
        pdbStruct.XYZ(i,1) = str2double(tline(31:38));
        pdbStruct.XYZ(i,2) = str2double(tline(39:46));
        pdbStruct.XYZ(i,3) = str2double(tline(47:54));
    end
    
    tline = fgetl(fid);
end

pdbStruct.chainSerNo(i+1:bufferSize) = [];
pdbStruct.chainID(i+1:bufferSize) = [];
pdbStruct.resSeq(i+1:bufferSize) = [];
pdbStruct.resName(i+1:bufferSize) = [];
pdbStruct.AtomSerNo(i+1:bufferSize) = [];
pdbStruct.AtomName(i+1:bufferSize) = [];
pdbStruct.XYZ(i+1:bufferSize,:) = [];

if(sum(isnan(pdbStruct.chainSerNo) | isinf(pdbStruct.chainSerNo)) > 0)
    error('Some chain serial numbers are not finite.');
end
if(sum(arrayfun(@(x)(strcmp(x,'')), pdbStruct.chainID)) > 0)
    error('Some chain identifiers are empty.');
end
if(sum(isnan(pdbStruct.resSeq) | isinf(pdbStruct.resSeq)) > 0)
    error('Some residue sequence numbers are not finite.');
end
if(sum(arrayfun(@(x)(strcmp(x,'')), pdbStruct.resName)) > 0)
    error('Some residue names are empty.');
end
if(sum(isnan(pdbStruct.AtomSerNo) | isinf(pdbStruct.AtomSerNo)) > 0)
    error('Some atom serial numbers are not finite.');
end
if(sum(arrayfun(@(x)(strcmp(x,'')), pdbStruct.AtomName)) > 0)
    error('Some atom names are empty.');
end
if(sum(sum(isnan(pdbStruct.XYZ) | isinf(pdbStruct.XYZ))) > 0)
    error('Some coordinates are not finite.');
end

fclose(fid);
end