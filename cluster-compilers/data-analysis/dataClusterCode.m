%% preliminary code

% reseed the random number generator; otherwise same random
% numbers are generated each time MATLAB restarts
rng shuffle 

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

% this won't overwrite anyway, but this is just a safety measure
% folder named "runName" is moved to "rundirpath_rel"
if exist(new_runDir,'dir') == 0 
    mkdir(new_runDir);        
end

% define fit routine
if any(strcmpi(fit_opt,{'bleach','bleaching','includeBleaching'}))
    fit_opt = 'bleach';
elseif any(strcmpi(fit_opt,{'noBleach','noBleaching','weakBleach',...
        'weakBleaching'}))
    fit_opt = 'nobleach';
else
    error('Unknown fitting routine specified.')
end

% tag for ics guess option
if use_ics_guess % use ics guess
    guess_tag = '--use-guess';
elseif ~use_ics_guess % no ics guess
    guess_tag = '--no-guess';
end
%

tauMax = max(tauVector);%%
runName = strcat('tauMax-',num2str(tauMax),'--kSqMax-',num2str(kSqMax),...
    '--',fit_opt,'-opt',guess_tag);
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
runName = new_runName;

movefile(new_jobDir,new_runDir); % move from "queued-jobs" folder
runDir = [new_runDir,filesep,runName]; % replace runDir definition

%% load data

load(loadPath)

analysisDir = [runDir,filesep,'analysis'];
mkdir(analysisDir)

% for convention
J = loadedMovie;
clear loadedMovie

%% fit intensity trace

[p,ci,bleach_fig] = bleachFit(J)
k_p_fit = p(2); % fitted value for bleaching rate

% save intensity trace with bleach profile fit
filename = [runDir filesep 'analysis' filesep runName '_intensity_trace.fig']; % save figure .fig
saveas(bleach_fig,filename)
filename = [runDir filesep 'analysis' filesep runName '_intensity_trace.pdf']; % save figure .pdf
saveas(bleach_fig,filename)

% store fit details in struct
bleach_fit_info.opt_params = p;
bleach_fit_info.ci = ci;

%% compute kICS autocorrelation

if redo_kics
    kSqVector = getKSqVector(J);
    [kSqVectorSubset,kSqSubsetInd] = getKSqVector(J,'kSqMin',kSqMin,'kSqMax',kSqMax);
    
    r_k_circ_uncut = zeros(length(kSqVector),length(tauVector));
    r_k_circ = zeros(length(kSqVectorSubset),length(tauVector));
    
    %
    tic
    
    r_k_norm = kICS3(J-repmat(mean(J,3),[1,1,size(J,3)])); % kICS autocorrelation function (ACF)
    if do_interp
        n_theta_arr = n_theta*ones(1,length(kSqVectorSubset));
        
        r_k_tau = r_k_norm(:,:,tauVector+1);
        parfor tau_i = 1:length(tauVector)
            r_k_circ(:,tau_i) = ellipticInterp(r_k_tau(:,:,tau_i),kSqVectorSubset,n_theta_arr);
        end
        delete(gcp)
    else
        r_k_circ_uncut = circular(r_k_norm(:,:,tauVector+1));
        r_k_circ = r_k_circ_uncut(kSqSubsetInd,:);
    end
    r_k_abs = abs(r_k_circ); % get rid of complex values
    
    toc
    %
    
    % if any(strcmpi(normByLag,{'none','noNorm',''})) % some normalization for when the kICS AC is not normalized, otherwise the fit is unreasonable
    %     max_value = max(max(r_k_abs));
    %     r_k_abs = r_k_abs/max_value;
    % end
    
    cd(runDir)
    filename = [runDir filesep 'analysis' filesep 'kICS_Data.mat'];
    save(filename,'kSqVector','r_k_abs','r_k_norm'); % save computed kICS ACF
end

%% plot kICS data

% for consistency 
kICSCorrSubset = r_k_abs; % corresponding subset of kICS ACF

ksq2plot = linspace(kSqVectorSubset(1),kSqVectorSubset(end),nPtsFitPlot); % |k|^2 for plotting best fit/theory curves

figure()
hold on
box on

% loop and plot over fixed time lag
color = lines(length(plotTauLags));
plotLegend = cell(1,length(plotTauLags));
h_sim_data = zeros(1,length(plotTauLags));
for tauInd = 1:length(plotTauLags) 
    % plot heuristic kICS ACF
    h_sim_data(tauInd) = plot(kSqVectorSubset,kICSCorrSubset(:,tauInd),...
        '.','markersize',16,'Color',color(tauInd,:)); 
    plotLegend{tauInd} = ['$\tau = ' num2str(plotTauLags(tauInd)) '$'];
end
% labeling
xlabel('$|\mathbf{k}|^2$ (pixels$^{-2}$)','interpreter','latex','fontsize',14)
ylabel('$\phi(|\mathbf{k}|^2,\tau)$','interpreter','latex','fontsize',14)
legend(h_sim_data,plotLegend,'fontsize',12,'interpreter','latex')
xlim([kSqVectorSubset(1) kSqVectorSubset(end)])
ylims = get(gca,'ylim');
ylim([0 ylims(2)])
tightfig(gcf) % no white-space (3rd party package; works for release 2015a)

% save plots without fits
filename = [analysisDir filesep runName '_nofit.fig']; % save figure .fig
saveas(gcf,filename)
filename = [analysisDir filesep runName '_nofit.pdf']; % save figure .pdf
saveas(gcf,filename)

%% ICS for guesses

if use_ics_guess
    ics_figpath = [runDir filesep 'analysis' filesep 'ics_figs'];
    ics_run = ICSCompiler(J,xi_lags,eta_lags,'subsets',n_ics_subs,...
        'saveFigs',ics_figpath);
    
    % save ICS info
    filename = [runDir filesep 'analysis' filesep 'ics_info.mat'];
    save(filename,'ics_run')
    
    % compute tau = 0 lag in kICS
    n_theta_arr = n_theta*ones(1,length(kSqVector));
    
    r_k_tau_0 = kICS3(J-repmat(mean(J,3),[1,1,size(J,3)]),'normByLag','none'); % kICS autocorrelation unnormalized
    r_k_circ_0 = ellipticInterp(r_k_tau_0(:,:,1),kSqVector,n_theta_arr);
    %
    
    % get param guesses for kICS fitting
    kics_figpath = [runDir filesep 'analysis'];
    [params_guess,lb,ub] = getkICSGuess(J,ics_run,p,kSqVector,r_k_circ_0,...
        ksq_noise_lb,'saveFigs',kics_figpath);
    %
    
    % save kICS guesses
    filename = [runDir filesep 'analysis' filesep 'kics_guesses.mat'];
    save(filename,'params_guess','lb','ub')
    %
end

%% kICS fitting

T = length(session_info.tRange);

% function handles of LSF error and fit functions
if strcmp(fit_opt,'bleach')
    % with bleaching
    err = @(params) timeIntkICSBleachFit(params,kSqVectorSubset,tauVector,...
        k_p_fit,T,'err',kICSCorrSubset);
    fit_fun = @(params,ksq,tau) timeIntkICSBleachFit(params,ksq,tau,...
        k_p_fit,T);
elseif strcmp(fit_opt,'nobleach')
    % without bleaching
    err = @(params) timeIntkICSFit(params,kSqVectorSubset,tauVector,...
        'err',kICSCorrSubset);
    fit_fun = @(params,ksq,tau) timeIntkICSFit(params,ksq,tau);
else
    error('Unknown fitting routine specified.')
end
%

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
%
% set initial guess to best global fit to get immediate
% convergence and correct Hessian from "fminunc.m"
x0 = opt_params; 

[x,~,~,~,~,hessian] = fminunc(err,x0);
disp(['params from "fminunc.m":',num2str(x)])
disp('variance on fit parameters:')
disp(num2str(inv(hessian)))
%

% store fit details in struct
kics_fit_info.opt_params = opt_params;
kics_fit_info.err_min = err_min;
kics_fit_info.hessian = hessian;

%% plot best fit curves

for tauInd = 1:length(plotTauLags) % loop and plot over fixed time lag
    plot(ksq2plot,fit_fun(opt_params,ksq2plot,plotTauLags(tauInd)),...
        'Color',color(tauInd,:)) % plot best fit kICS ACF
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
