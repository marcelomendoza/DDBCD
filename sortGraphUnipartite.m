function [A,Z,eta,perm,mu]=sortGraphUnipartite(A,Z,eta,mu)
% Function to sort graph accourding to clustering structure derived from
% the bi-partite IRM model. The clusters are sorted according to the cluster sizes.
%
% Usage:
%   [A,Z,eta,perm]=sortGraphUnipartite(A,Z,eta,mu)
%
% Input:
%    A      Adjacency matrix or multigraph
%    Z      Clustering assignment matrix
%    eta    matrix of between cluster relations
%    mu     (optional input)
%
% Output: 
%    A      sorted Adjacency matrix or multigraph
%    Z      sorted clustering assignment matrix
%    eta    sorted matrix of between cluster relations
%    perm   the permutation of the indices
%    mu     sorted cluster probabilityes (used for finite generative model)
%
% Written by Morten Mørup

% Sort Z and eta according to cluster size
sumZ=sum(Z,2);
[val,ind]=sort(sumZ,'descend');
Z=Z(ind,:);
if nargin>3
    mu=mu(ind);
else
    mu=[];
end
eta=eta(ind,:);
eta=eta(:,ind);

% Sort A and Z according to cluster assignment
noc=size(Z,1);
[val,perm]=sort((2.^(0:noc-1))*Z,'ascend');
Z=Z(:,perm);
if iscell(A)
    for n=1:length(A)
        A{n}=A{n}(perm,:);
        A{n}=A{n}(:,perm);
    end
else
    A=A(perm,:);
    A=A(:,perm);
end





