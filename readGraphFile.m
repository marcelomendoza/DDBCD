function [A]=readGraphFile(foo,IDMap)
% input: foo: content file (cites); IDMap: key -> id; (key (counter),id
% (fileID))
% output:  A: J x J adjacency matrix (double dtype)
J=length(IDMap);
A=sparse(J,J);
fid=fopen(foo);
line=fgetl(fid);
while ischar(line)
    tokens=split(line,' ');
    a=tokens{1};
    b=tokens{2};
    i=IDMap(a);
    j=IDMap(b);
    A(i,j)=1;
    A(j,i)=1;
    line=fgetl(fid);
end
fclose(fid);