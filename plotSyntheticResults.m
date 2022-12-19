function plotSyntheticResults_extended(A,West,Zestimated,G,eta_sorted)
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
disp(['AUC=' num2str(round(calcAUC(West,A)*100)/100)]);
