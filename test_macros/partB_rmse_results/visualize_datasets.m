% Visualizes synthetic data sets without missing values. Real centroids
% are illustrated as red circles.
clear, clc, close all;
addpath('../../datasets');
addpath('../../toolbox/general');
addpath('../../toolbox/clustering');
load('real_centers');
datasets = {'S1'; 'S2'; 'S3'; 'S4'; 'S5D2'; 'S2D2'; 'O200'; 'O2000'};
for i = 1:length(datasets)
    load(datasets{i});
    X = normalizedata(X);
    C = realCenters{i};
    figure;
    hold on;
    grid on;
    xlim([-1,1]);
    ylim([-1,1]);
    xticks(-1:0.5:1);
    yticks(-1:0.5:1);
    title(sprintf('%s',datasets{i}));
    scatter(X(:,1),X(:,2),'ko');
    scatter(C(:,1),C(:,2),'ro','filled');
end




