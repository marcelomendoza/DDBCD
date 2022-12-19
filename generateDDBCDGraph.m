function [A,Z,eta,gap,A_sorted,Z_sorted,eta_sorted,noc]=generateDDBCDGraph(J,D,nob,alpha,bp,bn,gap_prior,type,gap)
% Script to generate binary unipartite graphs according to the DD Bayesian Community Detection
% (BCD) generative model
%
% Usage:
%   [A,Z,eta,gap,A_sorted,Z_sorted,eta_sorted,noc]=generateDDBCDGraph(J,D,nob,alpha,bp,bn,gap_par,type)
% 
% Input:
%   J           size of graph 
%   D           distance matrix (square, JxJ)
%   nob         number of blobs induces by the attribute (used to compute a)
%   alpha       parameter for the Chinese Restaurant Process (CRP) 
%   bp          1x2 vector (default [10 1]) where 
%               bp(1): indicate within community prior link count
%               bp(2): indicate between community prior link count
%   bn          1x2 vector (default [1 10]) where 
%               bn(1): indicate within community prior non-link count
%               bn(2): indicate between community prior non-link count
%   gap_prior   prior distribution parameters for the gap parameter (default [1 1])
%   type        'UnDirected' (default) or 'Directed'
%
% Output:
%   A           Generated graph that has not been sorted, i.e. A=Bernoulli(Z'*eta*Z)
%   Z           noc x J generated assignment matrix, i.e. Z ~ Discrete(mu)
%   eta         noc x noc generated group relations, i.e. eta_{lm} ~ Beta(Bp,Bn)
%   gap         gap parameter
%   A_sorted    sorted J x J adjacency matrix of the generated graph,
%   Z_sorted    sorted noc x J generated assignment matrix
%   eta_sorted  sorted noc x noc generated group relations
%   noc         number of connected components in C

if nargin<8
    type='UnDirected';
end
if nargin<7
    gap_prior=[1 10];
end
if nargin<6
    bn=[1 10];    
end
if nargin<5
    bp=[10 1];
end
if nargin<4
    alpha=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate P (exponential decay)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=1./(nob+1);
S=exp(-D./a);
s=size(S);
idx=1:s(1)+1:s(1)*s(2);
S(idx)=S(idx)*alpha;
P=zeros(J,J);
for i=1:J
    row=S(i,:)./sum(S(i,:));
    P(i,:)=row;
end

%%%%%%%%%%%%%%%%%%%%%%%
% Generate assignments%
%%%%%%%%%%%%%%%%%%%%%%%
% We allocate data points conditioned on distance measurements

C=zeros(J,J);
for i=1:J % arrivals in order
   p=P(i,:);
   pp=cumsum(p);
   ind=find(rand<pp,1,'first'); % assignment according to pp
   C(i,ind)=1;      
end

%%%%%%%%%%%%%%%%%%%%
%Generate clusters %
%%%%%%%%%%%%%%%%%%%%
% We induce Z from C
G=digraph(C);
bins=conncomp(G, 'Type','weak');
Z=zeros(1,J);
for i=1:J
    z=bins(1,i);
    Z(z,i)=1;
end
noc=size(Z,1); % number of connected components

% Parameter for the beta prior imposed on eta
% We will assume the links drawn within and between communities are drawn
% from the same distribution specified by bp(1), bn(1) and bp(2), bn(2) respectively

% We next draw cluster relations eta from the Beta distribution specified
% by a and b
% Draw within community densities
eta_diag=betarnd(bp(1)*ones(noc,1),bn(1)*ones(noc,1));
%gap=betarnd(gap_prior(1),gap_prior(2)); 


% Draw between community densities
eta=zeros(noc);
for i=1:noc
    for ii=i+1:noc
        p=betainc(gap*min([eta_diag(i), eta_diag(ii)]),bp(2),bn(2)); % as bp(2)=bn(2)=1, p=x_lm 
        p=p.*rand;
        eta(i,ii) = betainv(p,bp(2),bn(2)); % as bp(2)=bn(2)=1, eta(i,ii)=p
    end
end
if strcmp(type,'UnDirected')
    eta=eta+eta'+diag(eta_diag); % eta + t(eta) + diag(eta)
end

% We finally draw links A from the Bernoulli likelihood
% We are currently interested in drawing an undirected graph thus we draw half of the links from
% the upper and half of the links from the lower triangular part of the
% adjacency matrix and add these parts together.
A = (Z'*eta*Z)>rand(J);
if strcmp(type,'UnDirected')
    A=triu(A,1);
    A=A+A';
    A=sparse(A);
end

[~,~,~]=sortGraphUnipartite(A,Z,eta,ones(1,J));

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display generated graph %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
subplot(2,3,1);
imshow(D);
title('Distance matrix (D)','FontWeight','Bold');
subplot(2,3,2);
plot(G);
title('Assignment graph (C)','FontWeight','Bold');
subplot(2,3,3);
mySpyPlot(A,[],1000/J);
title('Generated Graph','FontWeight','Bold')
subplot(2,3,4);
[A_sorted,Z_sorted,eta_sorted,~]=sortGraphUnipartite(A,Z,eta);
mySpyPlot(A_sorted,1000/J,Z_sorted);
%title('Sorted Generated Graph','FontWeight','Bold')
subplot(2,3,5);
imagesc(Z_sorted); colormap(1-gray); axis tight; title('Generated Clustering Z','FontWeight','Bold')
subplot(2,3,6);
imagesc(eta_sorted); colormap(1-gray); axis tight; title('Model Parameters \eta','FontWeight','Bold')



