function [A]=readGraphFilewebkb(foo,IDMap)
% input: foo: content file (cites); IDMap: key -> id; (key (counter),id
% (fileID)), sep = ' '
% output:  A: J x J adjacency matrix (double dtype)
J=length(IDMap);
A=sparse(J,J);
fid=fopen(foo);
line=fgetl(fid);
while ischar(line)
    tokens=split(line,' ');
    a=string(tokens{1});
    b=string(tokens{2});
    if isKey(IDMap,a) && isKey(IDMap,b)
        i=IDMap(a);
        j=IDMap(b);
        A(i,j)=1;
        A(j,i)=1;
    end
    line=fgetl(fid);
end
fclose(fid);