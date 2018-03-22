function dnaTop = assignSeqRand(dnaTop)

rng(0);     % Randon number seed, make the result reproduceable
nBase = numel(dnaTop);

% Initialize the sequence information for each base
for i = 1:nBase
    dnaTop(i).seq = '';
end

% Randomly generate the sequences
for i = 1:nBase
    if(strcmp(dnaTop(i).seq,''))
        dnaTop(i).seq = randNucleotide();
        
        i_comp = dnaTop(i).across;
        if(i_comp >= 0)     % Complementary base exists
            if(strcmp(dnaTop(i_comp).seq,''))   % Nucleotide in the complementary base is not assigned
                dnaTop(i_comp).seq = wspair(dnaTop(i).seq);
            end
        end
    end
end

% Check
for i = 1:nBase
    i_comp = dnaTop(i).across;
    if(i_comp >= 0)
        assert(strcmp(wspair(dnaTop(i).seq),dnaTop(i_comp).seq));
    end
end

end


% Randomly generate a nucleotide
function nt = randNucleotide()

ntSet = {'A', 'T', 'G', 'C'};
i = ceil(rand()*4);
nt = ntSet{i};

end


% Get the complementary base
function y = wspair(x)

if('A'==x || 'a'==x)
    y = 'T';
elseif('G'==x || 'g'==x)
    y = 'C';
elseif('C'==x || 'c'==x)
    y = 'G';
elseif('T'==x || 't'==x)
    y = 'A';
else
    error('Illegal base.');
end

end