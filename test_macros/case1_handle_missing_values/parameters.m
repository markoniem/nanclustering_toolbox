function params = parameters()
% Description: 
% Return initial parameters for test macros. 
%
% Output:
%     params - Parameters of test macro
%
params.dataset = 'iris';
params.myfuns = {@EED, @ESD, @PDS, @AD, @ICkNNI, @kNNI};
params.probmiss = [0.05, 0.15, 0.30, 0.60];
params.repetitions = 10;

end

