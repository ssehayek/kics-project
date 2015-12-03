%% Preliminary Code

rng shuffle % reseed the random number generator; otherwise same numbers are generated each time MATLAB restarts
disp(maxNumCompThreads)

run('getAnalysisInput') % gets all parameters from "analysisInput.txt.mv"

addpath(genpath(codePath))

paramStr = strcat('N_',num2str(N_diff),'_T_',num2str(T),'_D_',num2str(diffusion),'_k_on_',num2str(k_on),...
    '_k_off_',num2str(k_off),'_k_p_',num2str(k_p),'_SNR_',num2str(snr));
dirName = [saveSimDir paramStr];

if exist(dirName,'dir') ~= 0 && createSim % prevents overwriting saved data (should be improved as not all parameters are considered)
    disp('Simulations already exist for these parameters'); createSim = 0;
elseif ~createSim && ~fitSim % exits if no routine is chosen
    disp('"createSim" and "fitSim" cannot both be 0.'); exit
elseif fitSim && ~createSim && ~exist(dirName,'dir') %
    disp('Existing simulations do not exist.'); exit % exit if simulations don't exist and createSim=0
elseif exist(dirName,'dir') == 0 && createSim % make new simulation directory if creating a simulation
    mkdir(dirName)
end

%% Create Simulations

if createSim
    J = dronpaSim(sz,T,w0,N_diff,diffusion,k_on,k_off,k_p,...
        prob_agg,mean_agg_num,std_agg_dist,num_filaments,prob_place,'snr',snr);
    
    filename = [dirName filesep 'movie.mat'];
    save(filename,'J'); % save simulation movie
end

%% Run kICS

if fitSim && ~createSim
    filename = [dirName filesep 'movie.mat'];
    load(filename)
end

r_k_norm = kICS3(J-repmat(mean(J,3),[1,1,T]),'normByLag',normByLag); % kICS correlation of data
[r_k_circ,kSqVector] = circular(r_k_norm,floor(sz/2)-2); % circular average over |k|^2
r_k_abs = abs(r_k_circ); % get rid of complex values
if any(strcmpi(normByLag,{'none','noNorm',''})) % some normalization for when the kICS AC is not normalized, otherwise the fit is unreasonable
    max_value = max(max(r_k_abs));
    r_k_abs = r_k_abs/max_value;
end
cd(runDir)
mkdir('Analysis')

filename = [runDir 'Analysis/kICS_Data.mat'];
save(filename,'kSqVector','r_k_abs','r_k_norm'); % save computed kICS ACF

filename = [runDir 'Analysis/' paramStr '_nofit.fig']; % save figure .fig
saveas(gcf,filename)
filename = [runDir 'Analysis/' paramStr '_nofit.pdf']; % save figure .pdf
saveas(gcf,filename)

%% kICS Corr Plot

if fitSim
    kSqMinIndex = find(kSqVector >= kSqMin,1,'first') % lowest index, i, which satisfies kSqVector(i) >= kSqMin
    kSqMaxIndex = find(kSqVector <= kSqMax,1,'last') % highest index, j, which satisfies kSqVector(i) <= kSqMax
    kSqVectorSubset = kSqVector(kSqMinIndex:kSqMaxIndex); % all values satisfying kSqMin <= kSqVector <= kSqMax
    kSqSubsetInd = kSqMinIndex:kSqMaxIndex; % all indices satisfying kSqMin <= kSqVector(i) <= kSqMax
    
    kICSCorrSubset = r_k_abs(kSqSubsetInd,:); % corresponding subset of kICS ACF
    if legacy
        err = @(params) kICSNormTauFitFluctNoise(params,kSqVectorSubset,tauVector,'normByLag',normByLag,'err',kICSCorrSubset); % function handle of LSF error
    else
        err = @(params) kICSFitBiasFluct(params,kSqVectorSubset,tauVector,T,'normByLag',normByLag,'err',kICSCorrSubset,'useGPU',useGPU); % function handle of LSF error
    end
    
    tic
    if parallel
        parpool
        
        opts = optimoptions(@fmincon,'Algorithm','interior-point');
        problem = createOptimProblem('fmincon','objective',...
            err,'x0',params_guess,'lb',lb,'ub',ub,'options',opts);
        ms = MultiStart('UseParallel',true,'Display','final'); % FIXED: multi start object, the function is using previous options for 2013a adaptability
        [opt_params,err_min] = run(ms,problem,startPts)
        
        delete(gcp)
    else
        opts = optimoptions(@fmincon,'Algorithm','interior-point');
        problem = createOptimProblem('fmincon','objective',...
            err,'x0',params_guess,'lb',lb,'ub',ub,'options',opts);
        gs = GlobalSearch; % global search object
        [opt_params,err_min] = run(gs,problem)
    end
    fitTime = toc; disp(['fitTime = ' num2str(fitTime)]);
    
    kSqFit = linspace(kSqVectorSubset(1),kSqVectorSubset(end),nPtsFitPlot); % |k|^2 for plotting best fit function
    
    figure()
    hold on
    box on
    
    color = lines(length(plotTauLags)); % lines(length(plotTauLags));
    plotLegend = cell(1,length(plotTauLags));
    h_sim_data = zeros(1,length(plotTauLags));
    for tauInd = 1:length(plotTauLags) % loop and plot over fixed time lag
        if legacy
            plot(kSqFit,kICSNormTauFitFluctNoise(opt_params,kSqFit,plotTauLags(tauInd),'normByLag',normByLag),'Color',color(tauInd,:)) % plot best fit kICS ACF
        else
            plot(kSqFit,kICSFitBiasFluct(opt_params,kSqFit,plotTauLags(tauInd),T,'normByLag',normByLag,'useGPU',useGPU),'Color',color(tauInd,:)) % plot best fit kICS ACF
        end
        h_sim_data(tauInd) = plot(kSqVectorSubset,kICSCorrSubset(:,plotTauLags(tauInd)+1),'.','markersize',16,'Color',color(tauInd,:)); % plot heuristic kICS ACF
        plotLegend{tauInd} = ['$\tau = ' num2str(plotTauLags(tauInd)) '$'];
    end
    
    % title('kICS Autocorrelation Function')
    xlabel('$|\mathbf{k}|^2$ (pixels$^{-2}$)','interpreter','latex','fontsize',14)
    ylabel('$\phi(|\mathbf{k}|^2,\tau)$','interpreter','latex','fontsize',14)
    legend(h_sim_data,plotLegend,'fontsize',12,'interpreter','latex')
    xlim([kSqMin kSqMax])
    ylims = get(gca,'ylim');
    ylim([0 ylims(2)])
    tightfig(gcf)
    
    filename = [runDir 'Analysis/' paramStr '.fig']; % save figure .fig
    saveas(gcf,filename)
    filename = [runDir 'Analysis/' paramStr '.pdf']; % save figure .pdf
    saveas(gcf,filename)
    
    filename = [runDir 'Analysis/opt_params.mat'];
    save(filename,'opt_params','err_min');
end