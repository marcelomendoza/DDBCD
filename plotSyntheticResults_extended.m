function plotSyntheticResults_extended(A,West,Zestimated,G,eta_sorted,L,NOC,eta_,gap_)
% function to visualize the analysis of synthetically generated IRM data
%
% Usage:
%    plotSyntheticResults(A,West,Zestimated)
%
% Input:
%   A            Generated adjancency matrix
%   West         Link-predicted values
%   Zestimated   Estimated assignment matrix
%

J=size(A,1);
figure;
subplot(2,2,1);
imagesc(eta_sorted); colormap(1-gray); axis tight; title('Model Parameters \eta','FontWeight','Bold')
subplot(2,2,2);
[A_sorted,Z_Est_sorted]=sortGraphUnipartite(A,Zestimated,rand(size(Zestimated,1)));
imagesc(A_sorted); colormap(1-gray); axis off; title('A (sorted)','FontWeight','Bold')
subplot(2,2,3);
imagesc(Z_Est_sorted); colormap(1-gray); axis off; title('Estimated (sorted) Matrix Z','FontWeight','Bold')
subplot(2,2,4);
plot(G);
title('Assignment graph (C)','FontWeight','Bold');
max=length(L);
x=zeros(1,max);
y=zeros(1,max);
for i=1:max
    x(i)=i;
    y(i)=5;
end
figure;
subplot(4,1,1);
plot(L);
title('Log likelihood','FontWeight','Bold');
subplot(4,1,2);
plot(x,NOC,x,y);
%50 for webkb
ylim([0,150]);
title('Number of clusters','FontWeight','Bold');
subplot(4,1,3);
plot(eta_);
ylim([0,1]);
title('\eta_{1,1}','FontWeight','Bold');
subplot(4,1,4);
plot(gap_);
ylim([0,1]);
title('\gamma','FontWeight','Bold');
