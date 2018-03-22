function dnaTop = modify_dnaTop(dnaTop)

% Flip the up/down neighbors
for i = 1:numel(dnaTop)
    tmp = dnaTop(i).up;
    dnaTop(i).up = dnaTop(i).down;
    dnaTop(i).down = tmp;
end

end