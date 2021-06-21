% Description: 
% Evaluate performance of distance computation algorithms with missing data. 
% Capabilities of available data, partial distance, expected squared euclidean, 
% and expected euclidean strategies can be evaluated. Results will be saved to file.  
%
% Some test data sets were selected from UCI repository 
% (https://archive.ics.uci.edu/ml/index.php) for testing purpose 
% and are available in 'benchmark_data' folder  
%
% Obtained results are saved to 'results.xlsx' file
%
% Reference results are available in the following article: 
% E. Eirola, G. Doquire, M. Verleysen, and A. Lendasse, "Distance estimation in
% numerical data sets with missing values", Information Sciences, vol. 240, pp.
% 115–128, 2013.
%
% Inputs (see current values from 'parameters.m' file):
%             dataset - Selected data set 
%              myfuns - Used distance estimation function(s) 
%                       Multiple algorithms can be selected
%                       Alternatives: 
%                       @AD - Available data 
%                       @PDS - Partial distance 
%                       @ESD - Expected squared Euclidean
%                       @EED - Expected Euclidean  
%                       @ICkNNI - Incomplete-Case kNN Imputation
%                       @kNNI - kNN Imputation
%            probmiss - Probabilities of missing values generated 
%                       to input data set 
%         repetitions - The number of repetitions of the algorithm  
%
% Output:
%            Average values and (standard deviations) of rmse errors
%            saved to .mat file
%
clear, clc, close all;
addpath('../../benchmark_data/UCI_data/');
addpath('../../toolbox/preprocess');
addpath('../../toolbox/distance_strategies');
params = parameters();
dataset = params.dataset;
myfuns = params.myfuns; 
probmiss = params.probmiss;
repetitions = params.repetitions;
load(dataset);
X = normalizedata(X,'zscore');
N = size(X,1);
distFull = pdist(X)';
rmse = cell(length(myfuns),length(probmiss));
fprintf('Data set: %s\n',dataset);
fprintf('Selected missing values: ');
for i = 1:length(probmiss) 
    if mod(i,length(probmiss)) ~= 0, fprintf('%d%%, ',100*probmiss(i));
    else, fprintf('%d%%\n',100*probmiss(i)); end
end
fprintf('Number of repetitions: %d\n',repetitions);
fprintf('Performing distance computation...\n');
rmse_errors = zeros(repetitions,length(probmiss),length(myfuns));
for j = 1:length(myfuns)
    fprintf('Algorithm: %d/%d\n',j,length(myfuns));
    for m = 1:length(probmiss)
        rng(1);
        err = zeros(repetitions,1);
        for i = 1:repetitions
            Xm = genmissdata(X,probmiss(m));
            M = sum(any(isnan(Xm),2));
            myfun = myfuns{j};
            distMiss = myfun(Xm);
            distMiss(isnan(distMiss)) = mean(distMiss,'omitnan');
            err(i) = rmsefunc(distFull,distMiss,M,N);
        end 
        rmse_errors(:,j,m) =  err;
        rmse{j,m} = sprintf('%.3f(%.3f)',mean(err,'omitnan'),std(err,'omitnan'));
    end
end
fprintf('Done!\n');
save(sprintf('rmse_%s',dataset),'rmse');
save(sprintf('errors_%s',dataset),'rmse_errors');
rmpath('../../benchmark_data/UCI_data/');
rmpath('../../toolbox/preprocess');
rmpath('../../toolbox/distance_strategies');

