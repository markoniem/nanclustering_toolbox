% Visualizes results obtained using key point selection in an
% initialization phase of clustering. Results images are saved to img
% folded.
clear, clc, close all;
load('new_datamatrices');
load('centers');
maxClusters = 20;
datasets = {'S1'; 'S2'; 'S3'; 'S4'; 'S5D2'; 'S2D2'; 'O200'; 'O2000'; 'ecoli'; 'iris'; 'seeds'};
correctResults = [15, 15, 15, 15, 5, 2, 5, 5, 3, 2, 3];
probmiss = [0.00; 0.05; 0.10; 0.20];
%
for i = 1:length(datasets)
    for j = 1:length(probmiss)
       X = new_datamatrices{i,j};
       C = centers{correctResults(i)-1,i,j};
       figure;
       hold on;
       title(sprintf('%s (missing %d%%)',datasets{i},100*probmiss(j)))
       scatter(X(:,1),X(:,2),'k.');
       scatter(C(:,1),C(:,2),'ro','filled');
       %saveas(gcf,sprintf('img/%s_missing_%d',datasets{i},100*probmiss(j)),'png');
    end
end

