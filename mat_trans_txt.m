%% ¼ÇÂ¼»ùÒò
clc
clear

g = load('C.mat');
C = g.C;
p = size(C,1);

fid = fopen('gene_test.txt','w');
for k = 1:p
   
    x = C{k,4};
    [n,c] = size(x);
    fprintf(fid,'Symbol');
    fprintf(fid,'\n');
    for i = 1:n
        for j = 1:c
            a = char(x{i,j});
            fprintf(fid,'%s',a);
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
end
fclose(fid);