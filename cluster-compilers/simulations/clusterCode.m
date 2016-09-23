%% preliminary code

% reseed the random number generator; otherwise same random
% numbers are generated each time MATLAB restarts
rng shuffle

run('getAnalysisInput') % gets all parameters from "analysisInput.txt.mv"

addpath(genpath(codePath))
% cd to main project folder (this is "BASEDIR"
% in preSubScript)
cd(['..',filesep,'..']) 

% organization for simulations

% this is meant to better organize files & speed up simulations
if num_filaments == 0 || prob_place == 0 % no filaments means no immobile particles, or aggregates
    prob_agg = 0; mean_agg_num = 0; std_agg_dist = 0;
    num_filaments = 0; prob_place = 0;
end
if prob_agg == 0 || mean_agg_num == 0 % no aggregates
    prob_agg = 0; mean_agg_num = 0;
end
if N_diff == 0 || D == 0 % no diffusers, or diffusion
    N_diff = 0; D = 0;
end
if k_off == 0 % no blinking
    k_on = 1;
end
%

% strings with parameters used for naming
paramStr1 = strcat('D_',num2str(D),'_k_on_',num2str(k_on),...
    '_k_off_',num2str(k_off),'_k_p_',num2str(k_p)); % intrinsic properties of molecules string
paramStr2 = strcat('prob_agg_',num2str(prob_agg),'_mean_agg_num_',...
    num2str(mean_agg_num),'_std_agg_dist_',num2str(std_agg_dist)); % properties of aggregates string
paramStr3 = strcat('num_fils_',num2str(num_filaments),'_prob_place_',num2str(prob_place)); % properties of filaments string
if n_sub_frames ~= 1
    paramStr4 = strcat('N_diff_',num2str(N_diff),'_w0_',num2str(w0),'_sz_',num2str(sz),'_T_',...
        num2str(T),'_nsub_',num2str(n_sub_frames),'_SNR_',num2str(snr)); % intrinsic/extrinsic properties of experiment string (filename)
else
    paramStr4 = strcat('N_diff_',num2str(N_diff),'_w0_',num2str(w0),'_sz_',num2str(sz),'_T_',...
        num2str(T),'_SNR_',num2str(snr));
end
%

simfilepath_rel = [paramStr1,filesep,paramStr2,filesep,paramStr3,filesep,paramStr4]; % relative path to simulation file

simfilepath_abs = [saveSimDir,filesep,simfilepath_rel,'_rep_',num2str(loadRep),'.mat']  % absolute path to simulation file.
% this path is only complete if loadRep is non-empty,
% o.w. the path is partial.
% note that loadRep is empty if createSim=0.
simdirpath_rel = [paramStr1,filesep,paramStr2,filesep,paramStr3]; % relative path to simulation directory
simdirpath_abs = [saveSimDir,filesep,simdirpath_rel]; % relative path to simulation directory

if createSim
    if exist(simdirpath_abs,'dir') == 0 % create simualation dir if non-existing
        mkdir(simdirpath_abs);
    end
    % if simulation with these params already exists, create
    % another one
    [simfilepath_abs,num_rep] = iterateFilename(simfilepath_abs);
    disp(['creating simulation repetition: ',num2str(num_rep),'.'])
elseif ~createSim && ~fitSim % exits if no routine is chosen
    disp('"createSim" and "fitSim" cannot both be 0.'); 
    exit
elseif ~createSim && fitSim && ~exist(simfilepath_abs,'file') % exit if simulations don't exist and createSim=0
    disp('Simulations do not exist for these parameters. Please change value of createSim to create this simulation.'); 
    exit 
end

% organization for run folder
if ~isempty(loadRep) % this condition should only be true if createSim=0
    num_rep = loadRep;
    disp(['loading repetition ',num2str(num_rep),' for analysis.'])
end

tauMax = max(tauVector);
runName = strcat('tauMax-',num2str(tauMax),'--kSqMax-',num2str(kSqMax));
if isempty(runTag) == 0 % append rep number to end of "runName"
    new_runName = [runName,'--',runTag,'--rep-',num2str(num_rep)];  % new run name (with tag)
else
    new_runName = [runName,'--rep-',num2str(num_rep)]; % new run name (no tag)
end 
new_runDir = [baseDir,filesep,'queued-jobs',filesep,new_runName]; % append rep number to end of "runDir"
movefile(runDir,new_runDir) % rename queued job folder to "new_runDir" (within "queued-jobs" folder)

% redefine variables 
runDir = new_runDir; 
runName = new_runName;

if isempty(runClass) == 0
    rundirpath_rel = [runClass,filesep,simfilepath_rel]; % new relative run path (with class)
else
    rundirpath_rel = simfilepath_rel; % new relative run path (no class)
end

if exist(rundirpath_rel,'dir') == 0 % this won't overwrite anyway, but this is just a safety measure
    mkdir(rundirpath_rel);          % folder named "runName" is moved to "rundirpath_rel"
end

if exist([rundirpath_rel,filesep,runName],'dir') ~= 0 % exits if run with same name already exists
    disp('Run already exists under this name. Please change value of runName.');
    [error_dir,~] = iterateFilename(['errors',filesep,'error']);
    mkdir(error_dir)
    movefile(runDir,error_dir); % moves run to "error" folder
    exit
end

movefile(runDir,rundirpath_rel); % rename run directory
cd([rundirpath_rel,filesep,runName])
runDir = pwd; % replace runDir definition

%% create simulations

if createSim
    [J,sim_info] = dronpaSimParallel(sz,T,w0,N_diff,D,k_on,k_off,k_p,...
        prob_agg,mean_agg_num,std_agg_dist,num_filaments,prob_place,...
        'snr',snr,'subFrames',n_sub_frames); % create simulation
    
    save(simfilepath_abs,'J','sim_info'); % save simulation movie
elseif fitSim && ~createSim % if file wasn't created, load it
    load(simfilepath_abs)
end
% corrections to true parameters when fitting for mean on fraction, k_on/K,
% and sum of rates, K=k_on+k_off
K = sim_info.k_on + sim_info.k_off;

sim_info.eta_p = sim_info.eta_p*sim_info.k_on/K;

sim_info.true_params(2) = sim_info.k_on/K;
sim_info.true_params(3) = K;
sim_info.true_params(6) = sim_info.eta_p;
%

%% fit intensity trace

mkdir('analysis')

[p,ci,bleach_fig] = bleachFit(J)
k_p_fit = p(2); % fitted value for bleaching rate

% save intensity trace with bleach profile fit
filename = [runDir filesep 'analysis' filesep runName '_intensity_trace.fig']; % save figure .fig
saveas(bleach_fig,filename)
filename = [runDir filesep 'analysis' filesep runName '_intensity_trace.pdf']; % save figure .pdf
saveas(bleach_fig,filename)
close(gcf)

% store fit details in struct
bleach_fit_info.opt_params = p;
bleach_fit_info.ci = ci;

%% compute kICS autocorr

kSqVector = getKSqVector(J);
[kSqVectorSubset,kSqSubsetInd] = getKSqVector(J,'kSqMin',kSqMin,'kSqMax',kSqMax);

r_k_circ_uncut = zeros(length(kSqVector),length(tauVector));
r_k_circ = zeros(length(kSqVectorSubset),length(tauVector));

%
tic

r_k_norm = kICS3(J-repmat(mean(J,3),[1,1,size(J,3)]),'normByLag',normByLag); % kICS autocorrelation function (ACF)
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

%% ICS for guesses

if use_ics_guess && fitSim
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
    
    % check if guesses concur with true params
    param_names = {'diffusion','r','K','f_d','w0','eta_p'};
    for ii = 1:length(params_guess)
        if sim_info.true_params(ii) < lb(ii) || sim_info.true_params(ii) > ub(ii)
            warning(['True simulation param is not within guessed bounds',...
                ' for param: ''',param_names{ii},'''.'])
        end
    end
    %
end

%% kICS fitting

kICSCorrSubset = mean(abs(r_k_abs),3); % corresponding subset of kICS ACF

if fitSim
    % function handle of LSF error
    err = @(params) timeIntkICSFit(params,kSqVectorSubset,tauVector,...
        k_p_fit,T,'normByLag',normByLag,'err',kICSCorrSubset);
    
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
    disp(['true parameters: ',num2str(sim_info.true_params)])
    
    % get confidence intervals
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
end

%% plot simulation data

ksq2plot = linspace(kSqVectorSubset(1),kSqVectorSubset(end),nPtsFitPlot); % |k|^2 for plotting best fit/theory curves

figure()
hold on
box on

color = lines(length(plotTauLags));
plotLegend = cell(1,length(plotTauLags));
h_sim_data = zeros(1,length(plotTauLags));
for tauInd = 1:length(plotTauLags) % loop and plot over fixed time lag
    h_sim_data(tauInd) = plot(kSqVectorSubset,kICSCorrSubset(:,tauInd),...
        '.','markersize',16,'Color',color(tauInd,:)); % plot simulated kICS ACF
    plotLegend{tauInd} = ['$\tau = ' num2str(plotTauLags(tauInd)) '$'];
end
% labeling
xlabel('$|\mathbf{k}|^2$ (pixels$^{-2}$)','interpreter','latex','fontsize',14)
ylabel('$\phi(|\mathbf{k}|^2,\tau)$','interpreter','latex','fontsize',14)
legend(h_sim_data,plotLegend,'fontsize',12,'interpreter','latex')
xlim([kSqVectorSubset(1) kSqVectorSubset(end)])
ylims = get(gca,'ylim');
ylim([0 ylims(2)])
tightfig(gcf)

% save plots without fits
filename = [runDir filesep 'analysis' filesep runName '_nofit.fig']; % save figure .fig
saveas(gcf,filename)
filename = [runDir filesep 'analysis' filesep runName '_nofit.pdf']; % save figure .pdf
saveas(gcf,filename)

%% plot theoretical curves

% time-int theory
%
h_theory = zeros(1,length(plotTauLags));
for tauInd = 1:length(plotTauLags) % loop and plot over fixed time lag
    h_theory(tauInd) = plot(ksq2plot,timeIntkICSFit(sim_info.true_params,ksq2plot,plotTauLags(tauInd),k_p,T,'normByLag',normByLag),...
        '--','Color',color(tauInd,:)); % plot theoretical kICS ACF
end
tightfig(gcf)
%
% save plots with theory curves superimposed
filename = [runDir filesep 'analysis' filesep runName '_theory_time_int.fig']; % save figure .fig
saveas(gcf,filename)
filename = [runDir filesep 'analysis' filesep runName '_theory_time_int.pdf']; % save figure .pdf
saveas(gcf,filename)
%

% no time-int theory
%
delete(h_theory)
for tauInd = 1:length(plotTauLags) % loop and plot over fixed time lag
    h_theory(tauInd) = plot(ksq2plot,...
        kICSNormTauFitFluctNoiseBleachBlinkFrac(sim_info.true_params,...
        ksq2plot,plotTauLags(tauInd),k_p,T,'normByLag',normByLag),...
        '--','Color',color(tauInd,:)); % plot theoretical kICS ACF
end
tightfig(gcf)
% save plots with theory curves superimposed
filename = [runDir filesep 'analysis' filesep runName '_theory_no_time_int.fig']; % save figure .fig
saveas(gcf,filename)
filename = [runDir filesep 'analysis' filesep runName '_theory_no_time_int.pdf']; % save figure .pdf
saveas(gcf,filename)
%

%% plot best fit curves

if fitSim
    delete(h_theory) % delete theory curve handles
    
    for tauInd = 1:length(plotTauLags) % loop and plot over fixed time lag
        plot(ksq2plot,timeIntkICSFit(opt_params,ksq2plot,plotTauLags(tauInd),k_p_fit,T,'normByLag',normByLag),...
            'Color',color(tauInd,:)) % plot best fit kICS ACF
    end
    tightfig(gcf)
    
    % save plots with fits
    filename = [runDir filesep 'analysis' filesep runName '.fig']; % save figure .fig
    saveas(gcf,filename)
    filename = [runDir filesep 'analysis' filesep runName '.pdf']; % save figure .pdf
    saveas(gcf,filename)
    close(gcf)
    
    % save fit details
    filename = [runDir filesep 'analysis' filesep 'fit_info.mat'];
    save(filename,'bleach_fit_info','kics_fit_info');
    if output_mins
        save(filename,'kics_manymins','-append');
    end
end
