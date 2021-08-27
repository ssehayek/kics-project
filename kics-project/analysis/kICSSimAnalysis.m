% this code is for plotting and fitting a TICS blink/bleach autocorrelation
% including plots and fits for blink/bleach regimes

%% input variables

% time window
use_time_win = 1;
time_win = 200;

% logical for including windowing in fit model
win_fit = 0;

% tau range to fit
tauVector = 1:10;
% tau range to plot
plotTauLags = 1:5;

% min/max bounds on |k|^2
kSqMin = eps;
kSqMax = 1.8;

% min/max bounds for determining offset in tau=0 lag
ksq_min_noise = 10;
ksq_max_noise = 15;

% number of fits to compare
n_fits = 5;

% lower/upper bounds on fit parameters, [D,rho_on,K,p_D]
lb = eps*ones(1,4);
ub = [Inf,1,Inf,1];

% rois to consider (one per row; form [x0,y0,width,height])

roi = ...
    [1,1,64,64;
    65,65,64,64;
    1,65,64,64;
    65,1,64,64];

T = size(J,3);
nPtsFitPlot = 1e3;

%% main code

close all

mkdir('analysis')

for i = 1:size(roi,1)
    % roi image series
    J_roi = J(roi(i,2):roi(i,2)+roi(1,4)-1,roi(i,1):roi(i,1)+roi(i,3)-1,:);
    
    %%% compute kICS autocorr
    %
    tic
    
    % kICS autocorrelation function (ACF)
    if use_time_win
        time_win_varargin = {'timeWin',time_win};
    else
        time_win_varargin = {};
    end
    r_k = kICS(J_roi,'normByLag','none',time_win_varargin{:});
    
    r_k_0_sub = kICSSubNoise(r_k,ksq_min_noise,ksq_max_noise);
    
    r_k_0_circ = circular(r_k_0_sub(:,:,1));
    r_k_circ = circular(r_k);
    
    toc
    %
    %%% cut acf
    
    % get and cut |k|^2 vector
    [kSqVector,kSqInd] = getKSqVector(J_roi);
    [kSqVectorSubset,kSqSubsetInd] = getKSqVector(J_roi,'kSqMin',kSqMin,'kSqMax',kSqMax);
    
    % cut autocorrelation
    r_k_circ_cut = r_k_circ(kSqSubsetInd,tauVector+1);
    % cut normalization
    r_k_0_circ_cut = r_k_0_circ(kSqSubsetInd,1);
    
    % normalization
    r_k_norm = real(r_k_circ_cut./r_k_0_circ_cut);
    %
    %%% autocorrelation & fit plot
    %
    % fit function for entire autocorrelation
    if win_fit
        k_p = eps;
        
        err = @(params) kICSSlideWinFit(params,kSqVectorSubset,tauVector,k_p,T,time_win+1,...
            'err',r_k_norm,'symvars',{''});
        fit_fun = @(params,ksq,tau) kICSSlideWinFit(params,ksq,tau,k_p,T,time_win+1,'symvars',{''});
    else
        err = @(params) timeIntkICSFit(params,kSqVectorSubset,tauVector,...
            'err',r_k_norm,'symvars',{''});
        fit_fun = @(params,ksq,tau) timeIntkICSFit(params,ksq,tau,'symvars',{''});
    end
    
    opt_params = zeros(n_fits,4);
    err_min = zeros(n_fits,1);
    % fit n_fits times
    i_fit = 1;
    while i_fit <= n_fits
        try 
            params_guess = rand(1,4);
            opts = optimoptions(@fmincon,'Algorithm','interior-point');
            problem = createOptimProblem('fmincon','objective',...
                err,'x0',params_guess,'lb',lb,'ub',ub,'options',opts);
            gs = GlobalSearch; % global search object
            [opt_params(i_fit,:),err_min(i_fit)] = run(gs,problem);
            
            i_fit = i_fit + 1;
        end
    end
    [err_min,i_min] = min(err_min);
    opt_params = opt_params(i_min,:);
    %
    %%% plot simulation data
    %
    [~,plot_idx] = ismember(plotTauLags,tauVector);
    
    % |k|^2 for plotting best fit/theory curves
    ksq2plot = linspace(kSqVectorSubset(1),kSqVectorSubset(end),nPtsFitPlot);
    
    figure()
    hold on
    box on
    
    color = lines(length(plotTauLags));
    plotLegend = cell(1,length(plotTauLags));
    h_sim_data = zeros(1,length(plotTauLags));
    for tauInd = 1:length(plotTauLags) % loop and plot over fixed time lag
        h_sim_data(tauInd) = plot(kSqVectorSubset,r_k_norm(:,plot_idx(tauInd)),...
            '.','markersize',10,'Color',color(tauInd,:)); % plot simulated kICS ACF
        plotLegend{tauInd} = ['$\tau = ' num2str(plotTauLags(tauInd)) '$'];
    end
    % labeling
    xlabel('$|${\boldmath$k$}$|^2$ (pixels$^{-2}$)','interpreter','latex','fontsize',14)
    ylabel('$\tilde{\Phi}(|${\boldmath$k$}$|^2,\tau)$','interpreter','latex','fontsize',14)
    legend(h_sim_data,plotLegend,'fontsize',12,'interpreter','latex')
    xlim([kSqVectorSubset(1) kSqVectorSubset(end)])
    ylims = get(gca,'ylim');
    %
    %%% plot best fit curves
    %
    for tauInd = 1:length(plotTauLags)
        % plot best fit kICS ACF
        plot(ksq2plot,fit_fun(opt_params,ksq2plot,plotTauLags(tauInd)),...
            'Color',color(tauInd,:),'LineWidth',1.2)
    end
    tightfig(gcf)
    %
    %%% saving
    %
    roi_x1 = num2str(roi(i,1));
    roi_x2 = num2str(roi(i,1)+roi(i,3)-1);
    roi_y1 = num2str(roi(i,2));
    roi_y2 = num2str(roi(i,2)+roi(i,4)-1);
    filename = ['ROI_X_',roi_x1,'_',roi_x2,'_Y_',roi_y1,'_',roi_y2];
    
    save(['analysis/',filename],'r_k_norm','kSqMin','kSqMax','opt_params',...
        'win_fit','time_win','tauVector','err_min')
    %
end