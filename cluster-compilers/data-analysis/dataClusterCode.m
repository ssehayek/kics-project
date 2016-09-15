%% preliminary code

rng shuffle % reseed the random number generator; otherwise same random
            % numbers are generated each time MATLAB restarts

run('getAnalysisInput') % gets all parameters from "analysisInput.txt.mv"

addpath(genpath(codePath))

% organization for run folder

[loaddir,loadfilename,~] = fileparts(loadPath);

C = strsplit(loaddir,filesep);
loadClass = C{end}; % name of last directory under which the file is stored

C = strsplit(loadfilename,'.');
loadName = C{1}; % name of file without delimiters

if isempty(runClass) == 0
    new_runDir = [baseDir,filesep,'runs',filesep,runClass,filesep,loadClass,filesep,loadName]; % new absolute run path (with class)
else
    new_runDir = [baseDir,filesep,'runs',filesep,loadClass,filesep,loadName]; % new absolute run path (no class)
end

if exist(new_runDir,'dir') == 0 % this won't overwrite anyway, but this is just a safety measure
    mkdir(new_runDir);          % folder named "runName" is moved to "rundirpath_rel"
end

tauMax = max(tauVector);%%
runName = strcat('tauMax-',num2str(tauMax),'--kSqMax-',num2str(kSqMax));
if isempty(runTag) == 0 % append rep number to end of "runName"
    new_runName = [runName,'--',runTag];  % new run name (with tag)
else
    new_runName = runName; % new run name (no tag)
end 
new_fullpath = [new_runDir,filesep,new_runName,'--rep-']; 
[~,num_rep] = iterateFilename(new_fullpath); 
new_runName = [new_runName,'--rep-',num2str(num_rep)]; % append rep number at end of new filename

new_jobDir = [baseDir,filesep,'queued-jobs',filesep,new_runName]; % append "runTag" to end of "runDir"
movefile(runDir,new_jobDir) % rename queued job folder to string in "new_runDir" (within "queued-jobs" folder)

% redefine variables 
runName = new_runName;%%

movefile(new_jobDir,new_runDir); % move from "queued-jobs" folder
runDir = [new_runDir,filesep,runName]; % replace runDir definition

%% load data

load(loadPath)

analysisDir = [runDir,filesep,'analysis'];
mkdir(analysisDir)

%% fit intensity trace

[p,ci,bleach_fig] = bleachFit(loadedMovie)
k_p_fit = p(2); % fitted value for bleaching rate

% save intensity trace with bleach profile fit
filename = [runDir filesep 'analysis' filesep runName '_intensity_trace.fig']; % save figure .fig
saveas(bleach_fig,filename)
filename = [runDir filesep 'analysis' filesep runName '_intensity_trace.pdf']; % save figure .pdf
saveas(bleach_fig,filename)

% store fit details in struct
bleach_fit_info.opt_params = p;
bleach_fit_info.ci = ci;

%% fitting

kSqMinIndex = find(kSqVector >= kSqMin,1,'first') % lowest index, i, which satisfies kSqVector(i) >= kSqMin
if strcmpi(kSqMax,'max') % highest index, j, which satisfies kSqVector(j) <= kSqMax
    kSqMaxIndex = length(kSqVector); % use entire |k|^2 vector
else
    kSqMaxIndex = find(kSqVector <= kSqMax,1,'last') 
end
kSqVectorSubset = kSqVector(kSqMinIndex:kSqMaxIndex); % all values satisfying kSqMin <= kSqVector <= kSqMax
kSqSubsetInd = kSqMinIndex:kSqMaxIndex; % all indices satisfying kSqMin <= kSqVector(i) <= kSqMax

kICSCorrSubset = r_k_abs(kSqSubsetInd,:); % corresponding subset of kICS ACF

T = length(session_info.tRange);
err = @(params) kICSNormTauFitFluctNoiseBleach(params,kSqVectorSubset,tauVector,...
    k_p_fit,T,'normByLag',normByLag,'err',kICSCorrSubset); % function handle of LSF error

tic
if parallel
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
    if output_mins
        [opt_params,err_min,~,~,kics_manymins] = run(ms,problem,startPts);
    else
        [opt_params,err_min] = run(ms,problem,startPts);
    end

    delete(gcp) % delete parallel pool object
else
    opts = optimoptions(@fmincon,'Algorithm','interior-point');
    problem = createOptimProblem('fmincon','objective',...
        err,'x0',params_guess,'lb',lb,'ub',ub,'options',opts);
    gs = GlobalSearch; % global search object
    [opt_params,err_min] = run(gs,problem);
end
fitTime = toc; disp(['fitTime = ' num2str(fitTime)]);

disp(['optimal parameters: ',num2str(opt_params)])
disp(['minimum objective function: ',num2str(err_min)])
%     disp(['true parameters: ',num2str(sim_info.true_params)])

% get confidence intervals
x0 = opt_params; % set initial guess to best global fit to get immediate
                 % convergence and correct Hessian from "fminunc.m"

[x,~,~,~,~,hessian] = fminunc(err,x0);
disp(['params from "fminunc.m":',num2str(x)])
disp('variance on fit parameters:')
disp(num2str(inv(hessian)))
%

% store fit details in struct
kics_fit_info.opt_params = opt_params;
kics_fit_info.err_min = err_min;
kics_fit_info.hessian = hessian;

%% plot kICS data

ksq2plot = linspace(kSqVectorSubset(1),kSqVectorSubset(end),nPtsFitPlot); % |k|^2 for plotting best fit/theory curves

figure()
hold on
box on

color = lines(length(plotTauLags));
plotLegend = cell(1,length(plotTauLags));
h_sim_data = zeros(1,length(plotTauLags));
for tauInd = 1:length(plotTauLags) % loop and plot over fixed time lag
    h_sim_data(tauInd) = plot(kSqVectorSubset,kICSCorrSubset(:,plotTauLags(tauInd)+1),'.','markersize',16,'Color',color(tauInd,:)); % plot heuristic kICS ACF
    plotLegend{tauInd} = ['$\tau = ' num2str(plotTauLags(tauInd)) '$'];
end
% labeling
xlabel('$|\mathbf{k}|^2$ (pixels$^{-2}$)','interpreter','latex','fontsize',14)
ylabel('$\phi(|\mathbf{k}|^2,\tau)$','interpreter','latex','fontsize',14)
legend(h_sim_data,plotLegend,'fontsize',12,'interpreter','latex')
xlim([kSqMin kSqMax])
ylims = get(gca,'ylim');
ylim([0 ylims(2)])
tightfig(gcf) % no white-space (3rd party package; works for release 2015a)

% save plots without fits
filename = [analysisDir filesep runName '_nofit.fig']; % save figure .fig
saveas(gcf,filename)
filename = [analysisDir filesep runName '_nofit.pdf']; % save figure .pdf
saveas(gcf,filename)

%% plot best fit curves

for tauInd = 1:length(plotTauLags) % loop and plot over fixed time lag
    plot(ksq2plot,kICSNormTauFitFluctNoiseBleach(opt_params,ksq2plot,plotTauLags(tauInd),...
        k_p_fit,T,'normByLag',normByLag),'Color',color(tauInd,:)) % plot best fit kICS ACF
end
tightfig(gcf)

% save plots with fits
filename = [analysisDir filesep runName '.fig']; % save figure .fig
saveas(gcf,filename)
filename = [analysisDir filesep runName '.pdf']; % save figure .pdf
saveas(gcf,filename)

% save fit details
filename = [analysisDir filesep 'fit_info.mat'];
if output_mins
    save(filename,'bleach_fit_info','kics_fit_info','kics_manymins','-v7.3');
else
    save(filename,'bleach_fit_info','kics_fit_info');
end
