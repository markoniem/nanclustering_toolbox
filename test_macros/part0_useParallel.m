% Description:
%
% Performs weighted K-spatial-medians clustering using parallelled computation. 
% Validates results with replicated data set
%
clear, clc, close all;
rng(1); % for reproducibility purposes
addpath('../datasets');
addpath('../toolbox/general');
addpath('../toolbox/clustering');
addpath('../toolbox/cluster_validation');
params = parameters();
dataset = params.datasets;
probmiss = params.probmiss;
treatMissing = params.treatMissing;
treatMissingValidation = params.treatMissingValidation;
distance = params.distance;
pipelined = params.pipelined;
clusterInd = params.clusterInd;
useParallel = params.useParallel;
%
load(dataset{1,1});
[Xm, I] = genmissdata(X, probmiss);
Xm = Xm(I,:);
X = X(I,:);
weights = randi(5,size(X,1),1);
Xm_new = zeros(sum(weights,1),size(X,2));
for i = 1:length(weights)
    Xm_new = [Xm_new; repmat(Xm(i,:),weights(i),1)];
end
%
[centers, labels] = iterative_kcentroids(Xm, 'Weights', weights, 'TreatMissing', treatMissing, 'Distance', distance, ... 
                            'Pipelined', pipelined, 'UseParallel', useParallel, 'CloseParallel', false);
[indices_values, indices_names] = cluster_validation(Xm, centers, labels, 'Weights', weights, 'ClusterInd', clusterInd, ... 
                            'TreatMissing', treatMissingValidation, 'Distance', distance, 'UseParallel', useParallel);
[centers2, labels2] = iterative_kcentroids(Xm_new, 'TreatMissing', treatMissing, 'Distance', distance, ...
                            'Pipelined', pipelined, 'UseParallel', useParallel, 'CloseParallel', false);
[indices_values2, indices_names2] = cluster_validation(Xm_new, centers2, labels2, 'ClusterInd', clusterInd, ...
                             'TreatMissing', treatMissingValidation, 'Distance', distance, 'UseParallel', useParallel);
%                         
for i = 1:length(clusterInd)
    figure;
    hold on;
    subplot(2,1,1);
    plot(indices_values(i,:));
    title(char(indices_names{i}));
    subplot(2,1,2);
    plot(indices_values2(i,:));
end


function params = parameters()
% Description: 
% Returns initial parameters 
%
% Output:
%     params - Parameters of test macro
%
params.datasets = {'S4'};
params.probmiss = 0.10;
params.treatMissing = 'exp';
params.treatMissingValidation = 'exp';
params.distance = 'euc';
params.pipelined = false;
params.clusterInd = {@CalinskiHarabasz; @DaviesBouldin; @DaviesBouldin2; ... 
    @GenDunn; @kCE; @PBM; @RayTuri; @WB; @WemmertGancarski};
params.useParallel = true;

end



