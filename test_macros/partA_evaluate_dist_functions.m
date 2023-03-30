% Description: 
% Evaluate performance of distance computation algorithms with missing data. 
% Performance of available data, partial distance, expected squared euclidean, 
% expected euclidean, kNN-imputation, ICkNN-imputation, and matrix complementation 
% strategies are evaluated. Results are available in a xlsx file.  
%
% Some test data sets were selected from UCI repository 
% (https://archive.ics.uci.edu/ml/index.php) for testing purpose 
% and are available in 'dataset' folder  
%
% Reference results are available in the following article: 
% E. Eirola, G. Doquire, M. Verleysen, and A. Lendasse, "Distance estimation in
% numerical data sets with missing values", Information Sciences, vol. 240, pp.
% 115–128, 2013.
%
% Inputs: 
%       See parameters function
%
% Output:
%       Computed RMSE errors between real distances and estimated distances
%
clear, clc, close all;
addpath('../datasets');
addpath('../toolbox/general');
addpath('../toolbox/distance_computation');
%
params = parameters();
datasets = params.datasets;
myfuns = params.myfuns; 
probmiss = params.probmiss;
repetitions = params.repetitions;
%
% Set warnings off
warning off;
%
rmse = cell(length(myfuns),length(probmiss),length(datasets));
fprintf('Number of repetitions: %d\n',repetitions);
fprintf('Performing distance computation...\n');
%rmse_errors = zeros(repetitions,length(myfuns),length(probmiss));
for s = 1:length(datasets)
    fprintf('Data set: %s\n',datasets{s});
    load(datasets{s});
    X = normalizedata(X,'Method','zscore');
    for j = 1:length(myfuns)
        fprintf('Algorithm: %d/%d\n',j,length(myfuns));
        for m = 1:length(probmiss)
            fprintf('Percentage of missing values: %d %%\n',100*probmiss(m));
            rng(1);
            err = zeros(repetitions,1);
            for i = 1:repetitions
                [Xm,I] = genmissdata(X,probmiss(m));
                Xm = Xm(I,:);
                distFull = pdist(X(I,:))';
                N = size(Xm,1);
                M = sum(any(isnan(Xm),2));
                myfun = myfuns{j};
                [distMiss, success] = myfun(Xm);
                if ~success, err(i) = NaN; continue; end
                distMiss(isnan(distMiss)) = mean(distMiss,'omitnan');
                err(i) = rmsefunc(distFull,distMiss,M,N);
            end
            %rmse_errors(:,j,m) =  err;
            rmse{j,m,s} = sprintf('%.3f(%.3f)',mean(err,'omitnan'),std(err,'omitnan'));
        end
    end
    
end
fprintf('Done!\n');
restoredefaultpath;
%
% Save results to .xlsx file
%
distanceEstchar = cell(length(myfuns),1);
probMisschar = cell(length(probmiss),1);
for i = 1:length(myfuns), distanceEstchar{i} = sprintf('%s',char(myfuns{i,1})); end
for i = 1:length(probmiss), probMisschar{i} = sprintf('%s%%',string(100*probmiss(i))); end
for i = 1:size(rmse,3)    
      R = rmse(:,:,i);
      T = table(distanceEstchar,R);
      resultsTable = splitvars(T);
      resultsTable.Properties.VariableNames = [' ', probMisschar'];
      writetable(resultsTable,'partA_Results.xlsx', 'Sheet', datasets{i});
end

function params = parameters()
% Description: 
% Return initial parameters for test macro. 
%
% Output:
%     params - Parameters of test macro
%
params.datasets = {'iris'; 'ecoli'; 'breast_tissue'; 'glass'; 'wine'; 'parkinsons'; 'sonar'};
params.myfuns = {@ADS; @PDS; @ESD; @EED; @ICkNNI; @kNNI; @IST_MC; @ICkNNI_EXP; @kNNI_EXP};
params.probmiss = [0.05, 0.15, 0.30, 0.60];
params.repetitions = 250;

end


