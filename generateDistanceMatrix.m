function [D]=generateDistanceMatrix(J,noc)
% Script to generate a distance matrix according to a given number of blobs
%
% Usage:
%   [D]=generateDistanceMatrix(J,noc)
% 
% Input:
%   J           number of points
%   noc         number of blobs
%
% Output:
%   D           distance matrix (square, JxJ)
%
%%%%%%%%%%%%
%Generate X%
%%%%%%%%%%%%
delta=1./(noc+1);
params=zeros(2,noc);
for i=1:noc
    mu=i*delta;
    if mu<0.5
        alpha=(5*mu)/(1-mu);
        beta=5;
    elseif mu>0.5
        alpha=5;
        beta=(5-5*mu)/mu;
    else
        alpha=5;
        beta=5;
    end
    params(1,i)=alpha;
    params(2,i)=beta;
end
X=zeros(1,J);
pp=ones(1,noc);
pp=pp./noc;
p=cumsum(pp);
for i=1:J
    ind=find(rand<p,1,'first');
    alpha=params(1,ind);
    beta=params(2,ind);
    X(1,i)=betarnd(alpha,beta);
end
%%%%%%%%%%%%
%Generate D%
%%%%%%%%%%%%
D=zeros(J,J);
for i=1:J
    for j=1:J
        if i~=j
            D(i,j)=dist(X(1,i),X(1,j));
        end
    end
end
