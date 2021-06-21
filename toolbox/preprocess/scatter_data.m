function scatter_data(X)
% Description: 
% Visualize two-dimensional data.
%
% Function call:
%          scatter_data(X)
%
% Input:
%          X - Input data set 
%
figure;
hold on;
box on;
axis square;
title('Scattering result');
scatter(X(:,1),X(:,2),'bo','filled');

end

