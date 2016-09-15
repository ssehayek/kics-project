%% preliminary code

rng shuffle % reseed the random number generator; otherwise same numbers are generated each time MATLAB restarts

run('getAnalysisInput') % gets all parameters from "analysisInput.txt.mv"

addpath(genpath(codePath)) % add all directories under "codePath" to path variable

%% kICS logistics

[~,kSqVector] = circular(zeros(sz),floor(sz/2)-2); % |k|^2 values to examine

kSqMinIndex = find(kSqVector >= kSqMin,1,'first') % lowest index, i, which satisfies kSqVector(i) >= kSqMin
kSqMaxIndex = find(kSqVector <= kSqMax,1,'last') % highest index, j, which satisfies kSqVector(j) <= kSqMax
kSqVectorSubset = kSqVector(kSqMinIndex:kSqMaxIndex); % all values satisfying kSqMin <= kSqVector <= kSqMax
kSqSubsetInd = kSqMinIndex:kSqMaxIndex; % all indices satisfying kSqMin <= kSqVector(i) <= kSqMax

tauVector = 0:maxTau;

set_params = [D,k_on,k_off,f_d,w0,sigma]; % chosen params

% kICS function from set parameters
theory_pts = kICSNormTauFitFluctNoise(set_params,...
    kSqVectorSubset,tauVector,'normByLag',normByLag);

% define least squares error between fit and set curves
err = @(fit_params) kICSNormTauFitFluctNoise(fit_params,...
    kSqVectorSubset,tauVector,'normByLag',normByLag,'err',theory_pts);

%% fitting

tic
parpool % start parallel pool

opts = optimoptions(@fmincon,'Algorithm','interior-point');
problem = createOptimProblem('fmincon','objective',...
    err,'x0',params_guess,'lb',lb,'ub',ub,'options',opts);
ms = MultiStart('UseParallel',true,'Display','final'); 
ms.TolX = tolX; ms.TolFun = tolFun;

% scatter "startPts" many points (parallel) in parameter space and
% wait for local convergence (if possible) of each point.
% "opt_params" are the parameters yielding lowest value in
% objective function, "err_min"
[opt_params,err_min,~,~,manymins] = run(ms,problem,startPts);
opt_params
err_min

delete(gcp) % delete parallel pool object
fitTime = toc; disp(['fitTime = ' num2str(fitTime)]);

%% save results

mkdir([runDir,filesep,'analysis'])

filename = [runDir,filesep,'analysis',filesep,'manymins.mat'];
save(filename,'opt_params','err_min','manymins')
