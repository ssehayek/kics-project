%% preliminary code

% reseed the random number generator; otherwise same random
% numbers are generated each time MATLAB restarts
seed = rng('shuffle');

run('analysisInput.m') % gets all parameters from "analysisInput.m"

%% create save path for simulation

% this is meant to better organize files & speed up simulations
%
if (num_filaments == 0 || prob_place == 0) && (use_mask == 0)
    % no filaments means no immobile particles, or aggregates
    num_filaments = 0; prob_place = 0;
    mean_agg_num = 0; std_agg_dist = 0;
end
if use_mask
    num_filaments = 0; prob_place = 0;
end
if mean_agg_num == 0
    % no aggregates
    std_agg_dist = 0;
end
if N_diff == 0
    % no diffusers
    D = 0; frac_diff = 0;
end
if k_off == 0
    % no blinking
    k_on = 1;
end
%

dir_1 = ['D_',num2str(D),'_kon_',num2str(k_on),'_koff_',...
    num2str(k_off),'_kp_',num2str(k_p)];

if use_mask
    dir_2 = ['N_imm_',num2str(N_imm),'_mean_agg_num_',num2str(mean_agg_num),...
        '_std_agg_dist_',num2str(std_agg_dist),'_mask'];
else
    dir_2 = ['n_fils_',num2str(num_filaments),'_prob_place_',...
        num2str(prob_place),'_mean_agg_num_',num2str(mean_agg_num),...
        '_std_agg_dist_',num2str(std_agg_dist)];
end

if isempty(frac_diff)
    dir_3 = ['N_diff_',num2str(N_diff),'_w0_',...
        num2str(w0),'_sz_',num2str(sz),'_T_',num2str(T),'_nsub_',...
        num2str(nsub),'_noise_type_',num2str(noise_type)];
else
    dir_3 = ['frac_diff_',num2str(frac_diff),'_w0_',...
        num2str(w0),'_sz_',num2str(sz),'_T_',num2str(T),'_nsub_',...
        num2str(n_sub_frames),'_noise_type_',num2str(noise_type)];
end

if isempty(sim_tag)
    sim_tag = 'simulation';
end
    
save_simpath = iterateFilename([base_dir,filesep,'simulations',filesep...
    dir_1,filesep,dir_2,filesep,dir_3,filesep],'prepStr','_rep_');
save_simfile = [save_simpath,filesep,sim_tag,'.mat'];

if ~exist(save_simpath,'dir')
    mkdir(save_simpath)
end

disp(['creating simulation ',save_simfile,'.'])

%% create simulations

if createSim
    % struct for mask option
    mask_struct.use_mask = use_mask;
    mask_struct.mask_filepath = mask_filepath;
    mask_struct.N_imm = N_imm;
    
    [J,true_params] = kicsSim(sz,T,w0,N_diff,D,k_on,k_off,k_p,...
        mean_agg_num,std_agg_dist,num_filaments,prob_place,'parallel',sim_parallel,...
        'subFrames',n_sub_frames,'savepath',save_simfile,'noiseType',...
        noise_type,'kerVar',kernel_varargin,'laserVar',laser_varargin,...
        'noiseVar',noise_varargin,'fracDiff',frac_diff,'mask',mask_struct);
elseif fitSim % if simulation wasn't created, load it
    load(save_simfile)
end

%% fit intensity trace

[p,ci,bleach_fig] = bleachFit(J)
k_p_fit = p(2); % fitted value for bleaching rate

% if save_run
%     % save intensity trace with bleach profile fit
%     filename = [analysis_runpath filesep 'intensity_trace.fig']; % save figure .fig
%     saveas(bleach_fig,filename)
%     filename = [analysis_runpath filesep 'intensity_trace.pdf']; % save figure .pdf
%     saveas(bleach_fig,filename)
%     close(gcf)
% end

% store fit details in struct
bleach_fit_info.opt_params = p;
bleach_fit_info.ci = ci;

%% compute kICS autocorr

kSqVector = getKSqVector(J);
[kSqVectorSubset,kSqSubsetInd] = getKSqVector(J,'kSqMin',kSqMin,'kSqMax',kSqMax);

%
tic

% kICS autocorrelation function (ACF)
r_k = kICS(J,'normByLag','none');

r_k_circ = zeros(length(kSqVectorSubset),length(tauVector));

r_k_0_sub = kICSSubNoise(r_k,ksq_min_noise,ksq_max_noise);
if do_interp
    n_theta_arr = n_theta*ones(1,length(kSqVectorSubset));
    r_k_0_circ = ellipticInterp(r_k_0_sub,kSqVectorSubset,n_theta_arr);
    
    r_k_tau = r_k(:,:,tauVector+1);
    parfor tau_i = 1:length(tauVector)
        r_k_circ(:,tau_i) = ellipticInterp(r_k_tau(:,:,tau_i),kSqVectorSubset,n_theta_arr);
    end
    delete(gcp)
else
    r_k_0_circ_uncut = circular(r_k_0_sub(:,:,1));
    r_k_0_circ = r_k_0_circ_uncut(kSqSubsetInd,1);
    
    r_k_circ_uncut = circular(r_k(:,:,tauVector+1));
    r_k_circ = r_k_circ_uncut(kSqSubsetInd,:);
end
% normalization
r_k_norm = abs(r_k_circ)./abs(r_k_0_circ);
% r_k_norm = real(r_k_circ./r_k_0_circ(1,1));

toc
%

% if save_run
%     filename = [analysis_runpath filesep 'kics_data.mat'];
%     % save computed kICS ACF
%     save(filename,'kSqVector','r_k','r_k_norm');
% end

%% plot simulation data

% |k|^2 for plotting best fit/theory curves
ksq2plot = linspace(kSqVectorSubset(1),kSqVectorSubset(end),nPtsFitPlot);

figure()
hold on
box on

color = lines(length(plotTauLags));
plotLegend = cell(1,length(plotTauLags));
h_sim_data = zeros(1,length(plotTauLags));
for tauInd = 1:length(plotTauLags) % loop and plot over fixed time lag
    h_sim_data(tauInd) = plot(kSqVectorSubset,r_k_norm(:,tauInd),...
        '.','markersize',16,'Color',color(tauInd,:)); % plot simulated kICS ACF
    plotLegend{tauInd} = ['$\tau = ' num2str(plotTauLags(tauInd)) '$'];
end
% labeling
xlabel('$|\mathbf{k}|^2$ (pixels$^{-2}$)','interpreter','latex','fontsize',14)
ylabel('$\phi(|\mathbf{k}|^2,\tau)$','interpreter','latex','fontsize',14)
legend(h_sim_data,plotLegend,'fontsize',12,'interpreter','latex')
xlim([kSqVectorSubset(1) kSqVectorSubset(end)])
ylims = get(gca,'ylim');
% ylim([0 ylims(2)])
tightfig(gcf)

% if save_run
%     % save plots without fits
%     %
%     % save figure .fig
%     filename = [analysis_runpath filesep 'kics_acf.fig'];
%     saveas(gcf,filename)
%     % save figure .pdf
%     filename = [analysis_runpath filesep 'kics_acf.pdf'];
%     saveas(gcf,filename)
% end

%% plot theoretical curves

% % function handles of LSF error and fit functions
% err = @(params) kICSSlideWinImmBleachFit(params,kSqVectorSubset,tauVector,k_p,T,101,...
%     'err',r_k_norm,'symvars',{''});
% fit_fun = @(params,ksq,tau) kICSSlideWinImmBleachFit(params,ksq,tau,k_p,T,101,'symvars',{''});
err = @(params) timeIntkICSFit(params,kSqVectorSubset,tauVector,...
    'err',r_k_norm,'symvars',{''});
fit_fun = @(params,ksq,tau) timeIntkICSFit(params,ksq,tau,'symvars',{''});
% err = @(params) timeIntkICSBleachFit(params,kSqVectorSubset,tauVector,k_p_fit,T,...
%     'err',r_k_norm,'symvars',{''});
% fit_fun = @(params,ksq,tau) timeIntkICSBleachFit(params,ksq,tau,k_p,T,'symvars',{''});
% err = @(params) kICSSlideWinFullFit(params,kSqVectorSubset,tauVector,201,w0,...
%     'err',r_k_norm,'symvars',{''});
% fit_fun = @(params,ksq,tau) kICSSlideWinFullFit(params,ksq,tau,201,w0,'symvars',{''});
% err = @(params) kICSSlideWinFit(params,kSqVectorSubset,tauVector,k_p_fit,T,51,...
%     'err',r_k_norm,'symvars',{''});
% fit_fun = @(params,ksq,tau) kICSSlideWinFit(params,ksq,tau,k_p,T,51,'symvars',{''});

% err = @(params) kICSDiff3D(params,kSqVectorSubset,tauVector,...
%     'err',r_k_norm,'symvars',{''});
% fit_fun = @(params,ksq,tau) kICSDiff3D(params,ksq,tau,'symvars',{''});
% err = @(params) timeIntkICS3DFit(params,kSqVectorSubset,tauVector,...
%     'err',r_k_norm,'symvars',{''});
% fit_fun = @(params,ksq,tau) timeIntkICS3DFit(params,ksq,tau,'symvars',{''});
% err = @(params) kICSAnomDiff3D(params,kSqVectorSubset,tauVector,...
%     'err',r_k_norm,'symvars',{''});
% fit_fun = @(params,ksq,tau) kICSAnomDiff3D(params,ksq,tau,'symvars',{''});

% time-int theory
%
h_theory = zeros(1,length(plotTauLags));
for tauInd = 1:length(plotTauLags)
    h_theory(tauInd) = plot(ksq2plot,fit_fun(true_params,...
        ksq2plot,plotTauLags(tauInd)),...
        '--','Color',color(tauInd,:),'LineWidth',2);
end
tightfig(gcf)

% if save_run
%     % save plots with theory curves superimposed
%     filename = [analysis_runpath filesep 'kics_theory.fig']; % save figure .fig
%     saveas(gcf,filename)
%     filename = [analysis_runpath filesep 'kics_theory.pdf']; % save figure .pdf
%     saveas(gcf,filename)
% end

%% kICS fitting

if fitSim
    tic
    if fit_parallel
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
    disp(['true parameters: ',num2str(true_params)])
    
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

%% plot best fit curves

if fitSim
    delete(h_theory) % delete theory curve handles
    
    for tauInd = 1:length(plotTauLags)
        % plot best fit kICS ACF
        plot(ksq2plot,fit_fun(opt_params,ksq2plot,plotTauLags(tauInd)),...
            'Color',color(tauInd,:),'LineWidth',2)
    end
    tightfig(gcf)
    
    %     if save_run
    %         % save plots with fits
    %         %
    %         % save figure .fig
    %         filename = [analysis_runpath filesep 'kics_acf_fit.fig'];
    %         saveas(gcf,filename)
    %         % save figure .pdf
    %         filename = [analysis_runpath filesep 'kics_acf_fit.pdf'];
    %         saveas(gcf,filename)
    %         close(gcf)
    %         % save fit details
    %         filename = [analysis_runpath filesep 'fit_info.mat'];
    %         save(filename,'bleach_fit_info','kics_fit_info');
    %         if output_mins
    %             save(filename,'kics_manymins','-append','-v7.3');
    %         end
    %     end
end