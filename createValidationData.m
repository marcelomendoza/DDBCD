function [W,class]=createValidationData(A,pct,type)
% Create entries treated as missing in the graph(s) A
%
% Input:
% A     Cell array or sparse matrix of I x J adjacency matrices
% pct   percentage of links treated as missing - an equivalent number of
%       non-links will also be treated as missing in the generated W (default: pct=2.5)
% type  string indicating the following type of graphs
%         'UnDirected': Network is undirected
%         'Directed'  : Network is Directed
%         'Bipartite' : Network is Bipartite  (default)
%
% Output:
% W     Cell array or sparse matrix of I x J adjacency matrices indicating values to be used for validation
% class vector(s) indicating observed link and nonlinks in A stored in W,
%       i.e. class(k)=1 observation k in W is a link, class(k)=0,
%       observation k in W is a nonlink (for multi-graph this returns a cell array)

if nargin<3
    type='Bipartite';
end
if nargin<2
    pct=2.5;
end
pct=pct/100;

if iscell(A)
    [N1,N2]=size(A{1});
    for n=1:length(A)
        if pct==0
            W{n}=sparse(size(A{n},1),size(A{n},2));
        else
            if strcmp(type,'UnDirected')
                [I,J]=find(triu(A{n}));                
                density=length(I)/(N1*(N2-1)/2);
            else
                [I,J]=find(A{n});
                if strcmp(type,'Bipartite')
                    density=length(I)/(N1*N2);
                else
                    density=length(I)/(N1*(N2-1));
                end
            end
            M=length(I);
            ind=randperm(M);
            I1=I(ind(1:ceil(M*pct)));
            J1=J(ind(1:ceil(M*pct)));

            I2=ceil(size(A{n},1)*rand(ceil(M*pct*5),1));
            J2=ceil(size(A{n},2)*rand(ceil(M*pct*5),1));
            if strcmp(type,'UnDirected')
                ind=find(J2>I2);
                I2=I2(ind);
                J2=J2(ind);
            elseif strcmp(type,'Directed')
                ind=find(I2~=J2); % Make sure no diagonal elements treated as missing
                I2=I2(ind);
                J2=J2(ind);
            end
            I2=I2(1:min([ceil(M*pct/(1-density)),length(I2)]));
            J2=J2(1:min([ceil(M*pct/(1-density)),length(J2)]));        
            Wtmp=sparse(I2,J2,ones(1,length(I2)),N1,N2);
            Wtmp=Wtmp-Wtmp.*(A{n}>0);
            [I2,J2]=find(Wtmp);
            It=[I1; I2];
            Jt=[J1; J2];
            W{n}=sparse(It,Jt,ones(1,length(It)),size(A{n},1),size(A{n},2));
            [I,J]=find(W{n});
            W{n}=sparse(I,J,ones(1,length(I)),size(A{n},1),size(A{n},2));    
        end
        [~,~,class{n}]=find((A{n}+W{n}).*W{n});
        class{n}=class{n}-1;
    end
else
    [N1,N2]=size(A);
    if pct==0
        W=sparse(size(A,1),size(A,2));
        %class=[];
    else
        if strcmp(type,'UnDirected')
            [I,J]=find(triu(A));
            density=length(I)/(N1*(N2-1)/2);
        else
            [I,J]=find(A);
             if strcmp(type,'Bipartite')
                 density=length(I)/(N1*N2);
             else
                 density=length(I)/(N1*(N2-1));
             end
        end
        M=length(I);
        ind=randperm(M);
        I1=I(ind(1:ceil(M*pct)));
        J1=J(ind(1:ceil(M*pct)));

        I2=ceil(size(A,1)*rand(ceil(M*pct*5),1));
        J2=ceil(size(A,2)*rand(ceil(M*pct*5),1));
        if strcmp(type,'UnDirected')
            ind=find(J2>I2);
            I2=I2(ind);
            J2=J2(ind);
        elseif strcmp(type,'Directed')
            ind=find(I2~=J2); % Make sure no diagonal elements treated as missing
            I2=I2(ind);
            J2=J2(ind);
        end
        I2=I2(1:min([ceil(M*pct/(1-density)),length(I2)]));
        J2=J2(1:min([ceil(M*pct/(1-density)),length(J2)]));        
        Wtmp=sparse(I2,J2,ones(1,length(I2)),N1,N2);
        Wtmp=Wtmp-Wtmp.*(A>0);
        [I2,J2]=find(Wtmp);
        It=[I1; I2];
        Jt=[J1; J2];
        W=sparse(It,Jt,ones(1,length(It)),size(A,1),size(A,2));
        [I,J]=find(W);
        W=sparse(I,J,ones(1,length(I)),size(A,1),size(A,2));    
    end
    [~,~,class]=find((A+W).*W);
    class=class-1;
end
