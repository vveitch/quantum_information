n=4;
M = stabilizers(n);

for i=1:size(M,4)
    for j=1:size(M,3)
            dlmwrite('d16_stabs.txt', M(:,:,j,i), '-append', ...
                'roffset', 1, 'delimiter', ' ')
    end
end