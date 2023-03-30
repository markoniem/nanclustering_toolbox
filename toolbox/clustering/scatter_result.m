function scatter_result(X, C, L)
% Description: 
% Visualizes clustering result for 2D data.
%
% Input:
%          X - Input data set 
%          C - Cluster centroids
%          L - Cluster labels 
%
if (size(X,2) ~= 2)
    error('2D data is required');
end
k = size(C,1); 
figure;
hold on;
box on;
axis square;
colors = lines(k);
title('Clustering result');
for i = 1:k
    scatter(X(L==i,1),X(L==i,2),'CData',colors(i,:));
end
scatter(C(:,1),C(:,2),50,'ko','filled');

end

