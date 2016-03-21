%% Preliminary Code

rng shuffle % reseed the random number generator; otherwise same random 
            % numbers are generated each time MATLAB restarts

run('getAnalysisInput') % gets all parameters from "analysisInput.txt.mv"

addpath(genpath(codePath))
cd(['..',filesep,'..']) % cd to main project folder (this is "BASEDIR" 
                        % in preSubScript)

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
paramStr4 = strcat('N_diff_',num2str(N_diff),'_w0_',num2str(w0),'_sz_',num2str(sz),'_T_',...
    num2str(T),'_SNR_',num2str(snr)); % intrinsic/extrinsic properties of experiment string (filename)
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
    [simfilepath_abs,num_rep] = iterateFilename(simfilepath_abs); % if simulation with these
                                                                  % params already exists, create
                                                                  % another one
    disp(['creating simulation repetition: ',num2str(num_rep),'.'])
elseif ~createSim && ~fitSim % exits if no routine is chosen
    disp('"createSim" and "fitSim" cannot both be 0.'); 
elseif ~createSim && fitSim && ~exist(simfilepath_abs,'file') % exit if simulations don't exist and createSim=0
    disp('Simulations do not exist for these parameters. Please change value of createSim to create this simulation.'); 
    exit 
end

% organization for run folder
if ~isempty(loadRep) % this condition should only be true if createSim=0
    num_rep = loadRep;
    disp(['loading repetition ',num2str(num_rep),' for analysis.'])
end

new_runName = [runName,'--rep-',num2str(num_rep)]; % append rep number to end of "runName"
new_runDir = [baseDir,filesep,'queued-jobs',filesep,new_runName]; % append rep number to end of "runDir"
movefile(runDir,new_runDir) % rename folder with rep number appended

% redefine variables to have rep number appended
runDir = new_runDir; 
runName = new_runName;

if isempty(runClass) == 0
    rundirpath_rel = [runClass,filesep,simdirpath_rel]; % new relative run path (with class)
else
    rundirpath_rel = simdirpath_rel; % new relative run path (no class)
end

if exist(rundirpath_rel,'dir') == 0 % this won't overwrite anyway, but this is just a safety measure
    mkdir(rundirpath_rel);          % folder named "runName" is moved to "rundirpath_rel"
end

if exist([rundirpath_rel,filesep,runName],'dir') ~= 0 % exits if run with same name already exists
    disp('Run already exists under this name. Please change value of runName.');
    [error_dir,~] = iterateFilename('error');
    movefile(runDir,error_dir); % moves run to "error" folder
    exit
end

movefile(runDir,rundirpath_rel); % rename run directory
cd([rundirpath_rel,filesep,runName])
runDir = pwd; % replace runDir definition

%% Create Simulations

if createSim
    [J,sim_info] = dronpaSim(sz,T,w0,N_diff,D,k_on,k_off,k_p,...
        prob_agg,mean_agg_num,std_agg_dist,num_filaments,prob_place,'snr',snr); % create simulation
    
    save(simfilepath_abs,'J','sim_info'); % save simulation movie
end

%% Run kICS

if fitSim && ~createSim % if file wasn't created, load it
    load(simfilepath_abs)
end

r_k_norm = kICS3(J-repmat(mean(J,3),[1,1,T]),'normByLag',normByLag); % kICS correlation of data
[r_k_circ,kSqVector] = circular(r_k_norm); % circular average over |k|^2
r_k_abs = abs(r_k_circ); % get rid of complex values
if any(strcmpi(normByLag,{'none','noNorm',''})) % some normalization for when the kICS AC is not normalized, otherwise the fit is unreasonable
    max_value = max(max(r_k_abs));
    r_k_abs = r_k_abs/max_value;
end
cd(runDir)
mkdir('analysis')

filename = [runDir filesep 'analysis' filesep 'kICS_Data.mat'];
save(filename,'kSqVector','r_k_abs','r_k_norm'); % save computed kICS ACF

%% kICS Corr Plot

if fitSim
    kSqMinIndex = find(kSqVector >= kSqMin,1,'first') % lowest index, i, which satisfies kSqVector(i) >= kSqMin
    kSqMaxIndex = find(kSqVector <= kSqMax,1,'last') % highest index, j, which satisfies kSqVector(j) <= kSqMax
    kSqVectorSubset = kSqVector(kSqMinIndex:kSqMaxIndex); % all values satisfying kSqMin <= kSqVector <= kSqMax
    kSqSubsetInd = kSqMinIndex:kSqMaxIndex; % all indices satisfying kSqMin <= kSqVector(i) <= kSqMax
    
    kICSCorrSubset = r_k_abs(kSqSubsetInd,:); % corresponding subset of kICS ACF
    err = @(params) kICSNormTauFitFluctNoise(params,kSqVectorSubset,tauVector,'normByLag',normByLag,'err',kICSCorrSubset); % function handle of LSF error
    
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
        [opt_params,err_min,~,~,manymins] = run(ms,problem,startPts);
        disp(['optimal parameters: ',num2str(opt_params)])
        disp(['minimum objective function: ',num2str(err_min)])
        disp(['true parameters: ',num2str(sim_info.true_params)])

        delete(gcp) % delete parallel pool object
    else
        opts = optimoptions(@fmincon,'Algorithm','interior-point');
        problem = createOptimProblem('fmincon','objective',...
            err,'x0',params_guess,'lb',lb,'ub',ub,'options',opts);
        gs = GlobalSearch; % global search object
        [opt_params,err_min] = run(gs,problem);
    end
    fitTime = toc; disp(['fitTime = ' num2str(fitTime)]);
    
    kSqFit = linspace(kSqVectorSubset(1),kSqVectorSubset(end),nPtsFitPlot); % |k|^2 for plotting best fit function result
    
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
    tightfig(gcf)
    
    % save plots without fits
    filename = [runDir filesep 'analysis' filesep paramStr4 '_nofit.fig']; % save figure .fig
    saveas(gcf,filename)
    filename = [runDir filesep 'analysis' filesep paramStr4 '_nofit.pdf']; % save figure .pdf
    saveas(gcf,filename)
    
    for tauInd = 1:length(plotTauLags) % loop and plot over fixed time lag
        plot(kSqFit,kICSNormTauFitFluctNoise(opt_params,kSqFit,plotTauLags(tauInd),'normByLag',normByLag),'Color',color(tauInd,:)) % plot best fit kICS ACF
    end
    tightfig(gcf)
    
    % save plots with fits
    filename = [runDir filesep 'analysis' filesep paramStr4 '.fig']; % save figure .fig
    saveas(gcf,filename)
    filename = [runDir filesep 'analysis' filesep paramStr4 '.pdf']; % save figure .pdf
    saveas(gcf,filename)
    
    filename = [runDir filesep 'analysis' filesep 'fit_info.mat'];
    save(filename,'opt_params','err_min','manymins','-v7.3');
end	

if output_mins
	  save(filename,'opt_params','err_min','manymins','-v7.3');
else
	  save(filename,'opt_params','err_min');
end