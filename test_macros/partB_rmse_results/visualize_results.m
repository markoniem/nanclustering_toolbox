% Visualize clustering results for data sets with 20 % missing values on data. 
% Missing values are regenerated and clustering is performed 100 times 
% (with 100 replicates in each time). 
%
% Clustering is performed four different manners:
% 1) Clustering is initialized with real centroids and available data
% strategy is used for treating missing values. 
% 2) Clustering is not initialized 
% and available data strategy is used for distance computation in clustering.
% 3) Clustering is not initialized and expected Euclidean distance is used for 
% distance computation in clustering.
% 4) Clustering is not initialized and clustering results based on 
% expected Euclidean distances are pipelined to the clustering method with 
% available data strategy.
%
clear, clc, close all;
addpath('../../datasets');
load('real_centers');
load('realCentersArr');
load('centroids');
datasets = {'S1'; 'S2'; 'S3'; 'S4'; 'S5D2'; 'S2D2'; 'O200'; 'O2000'};
repetitions = 100;
titles = {'(ads: real cent. init)','(ads)','(eed)','(eed-ads)'};
titles2 = {'ads_real_cent_init','ads','eed','eed_ads'};
%
for i = 1:length(datasets)
    for k = 1:4
        X = [];
        C = realCenters{i};
        for j = 1:repetitions
            if k==1
                X = [X; realCentersArr{i,j}];
            else
                X = [X; centroids{i,j,k-1}];
            end
        end
        figure;
        hold on;
        grid on;
        xlim([-1,1]);
        ylim([-1,1]);
        xticks(-1:0.5:1);
        yticks(-1:0.5:1);
        title(sprintf('%s %s',datasets{i},titles{k}));
        scatter(X(:,1),X(:,2),'ko');
        scatter(C(:,1),C(:,2),'ro','filled');
        %saveas(gcf,sprintf('img/png/%s_%s',datasets{i},titles2{k}),'png');
        %saveas(gcf,sprintf('img/fig/%s_%s',datasets{i},titles2{k}),'fig');
    end
end

