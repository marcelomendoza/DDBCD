function [IDs,IDMap,D,dim]=readContentFilewebkb(foo,eta)
% input: foo: content file webkb (e.g., cornell.content), sep = '\t'
% output:  IDs: <id,ID,label> pairs
%          IDMap: M:i->id
%          D: J x J distance matrix (double dtype) 
fid=fopen(foo);
line=fgetl(fid);
i=1;
while ischar(line)
    tokens=split(line,'	');
    url=tokens(1);
    label=tokens(length(tokens));
    tuple={i,url,label};
    if i==1
        IDs={tuple};
    else
        IDs{end+1}=tuple;
    end
    dim=length(tokens)-2;
    X=zeros(1,dim);
    lim=dim+1;
    for j=2:lim
        entry=str2double(tokens{j});
        X(j-1)=entry;
    end
    if i==1
        Xs={X};
    else
        Xs{end+1}=X;
    end
    line=fgetl(fid);
    i=i+1;
end
fclose(fid);
J=length(Xs);

%{Cornell, Wisconsin, Texas, Washington, Citeseer, ...}
%cosine%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D=zeros(J,J);
for i=1:J
    X1=Xs{i}/sum((Xs{i}).^2).^0.5;
    for j=1:J
        X2=Xs{j}/sum((Xs{j}).^2).^0.5;
        if i~=j
            D(i,j)=(1-cosineSimilarity(X1,X2))^eta;
        end
    end
end
Y=max(D);
Norm=max(Y);
D=D./Norm;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{synth1,synth2}
%euclidean 1D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%D=zeros(J,J);
%for i=1:J
%    X1=Xs{i};
%    for j=1:J
%        X2=Xs{j};
%        if i~=j
%            D(i,j)=sum((X1-X2).^2).^0.5;
%        end
%    end
%end
%Y=max(D);
%Norm=max(Y);
%D=D./Norm;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IDMap=containers.Map;
for i=1:J
    val=IDs{i}{1}; %counter
    key=string(IDs{i}{2}); %fileID
    IDMap(key)=val;
end