% Description: 
% Data pre-processing for data set mapping
clear, clc, close all; 
addpath('../../benchmark_data/M_spheres');
addpath('../../toolbox/preprocess');
params = params();
datasets = params.datasets;
probmiss = params.probmiss;
distance = params.distance;
for i = 1:length(datasets)
    load(datasets{i});
    clear C;
    clear labels;
    fprintf('Data set: %s\n',datasets{i});
    datamatrices = cell(1,length(probmiss));
    for j = 1:length(probmiss)
        fprintf('Missing values: %.2f %% \n',probmiss(j));  
        X = normalizedata(X,'min-max');
        Xm = genmissdata(X,probmiss(j));   
        K = 5;
        Ximp = ICknnimpute(Xm,K,distance);
        opts = statset('Display','iter');
        [~, start, ~] = pca(Ximp);
        datamatrices{1,j} = mdscale(pdist(Ximp),2,'Start',start(:,1:2),'Options',opts); 
    end
    save(sprintf('M_Shpheres2D/%s.mat',char(datasets{i})),'datamatrices');
end

