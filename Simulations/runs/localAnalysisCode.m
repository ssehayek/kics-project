rng shuffle % reseed the random number generator; otherwise same random 
            % numbers are generated each time MATLAB restarts

%% Input Variable

%%%%%% Internal Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

codePath = 'C:\Users\SimonS\Dropbox (Personal)\Research\PhD\SOFI-Project'; % path where all codes are stored (recursive)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Simulation Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (refer to dronpaSim.m for variable definitions)

% main params
sz = 64; 
T = 50; 
n_sub_frames = 10;
w0 = 0.61*522/(2.8*178); 
N_diff = 500; 
D = 1; 
k_on = 0.2; 
k_off = 0.1; 
k_p = 0; 

% aggregate params
prob_agg = 0.5;
mean_agg_num = 5;
std_agg_dist = 0.1; 

% filament params 
num_filaments = 20; 
prob_place = 0.3; 

% varargin
snr = 8; % signal-to-noise ratio (enter 'none' for no noise)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% ICS Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

use_ics_guess = 1;
%%% relevant only if value is 1 %%%

% split data into 'n_ics_subs' many subsets  
n_ics_subs = 10;
% spatial lags to fit in ICS
xi_lags = -15:15;
eta_lags = xi_lags;
%

% lower bound of |k|^2 for which kICS tau = 0 lag (unnormalized) is noise
% dominated (should appear flat)
% CHOOSE THIS VALUE CAREFULLY
ksq_noise_lb = 15;

%%%%%%

%%%%%% kICS Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% choose whether to Fourier interpolate while circularly averaging
do_interp = 1;

% number of angles to sample for fixed |k|
n_theta = 1000;

%%%%%% kICS Fitting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

do_fit = 0; % run fit algorithm
do_theory = 1; % superimpose theory curves on simulation curves

normByLag = 0; % actual lag values i.e. tau=0 is 0th lag

% fitting 
tauVector = 1:5; % tau values to fit

%%% OLD FORMAT
% format [D,k_on,k_off,f_d,w0,eta_p]
% here "f_d" is the fraction of diffusing particles, and "eta'" is a term
% which is related to the noise
%
% params_guess = 100*rand(1,6); % guess params for fit
% lb = eps*ones(1,6); % lower bound for fit params
% ub = [10 1 1 1 10 Inf]; % upper bound
%
%%% NEW FORMAT
% format [D,r,K,f_d,w0,eta_p]
% here "K" is the sum of the blinking rates, and "r" is the fraction
% k_on/K; the other terms are defined as before ("eta_p" now has an extra
% factor of "r")
params_guess = [10*rand,rand,2*rand,rand,10*rand,1000*rand]; % guess params for fit
lb = eps*ones(1,6); % lower bound for fit params
ub = [10,1,2,1,10,Inf]; % upper bound

% plot kICS ACF
plotTauLags = tauVector; % tau values to plot
nPtsFitPlot = 1000; % number of |k|^2 pts to plot in best fit function
% min |k|^2 value to fit/plot
kSqMin = 0.01; 
% max |k|^2 value to fit/plot, put 'max' to fit entire range
kSqMax = 'max'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Simulations

addpath(genpath(codePath)) % add all codes to path variable (recursive)

J = zeros(sz,sz,T); 
% create simulation
[J,sim_info] = dronpaSim(sz,T,w0,N_diff,D,k_on,k_off,k_p,...
    prob_agg,mean_agg_num,std_agg_dist,num_filaments,prob_place,'snr',snr,...
    'subframes',n_sub_frames);
%

% corrections to true parameters when fitting for mean on fraction, k_on/K,
% and sum of rates, K=k_on+k_off
K = k_on + k_off;

sim_info.eta_p = sim_info.eta_p*k_on/K;

sim_info.true_params(2) = k_on/K;
sim_info.true_params(3) = K;
sim_info.true_params(6) = sim_info.eta_p;
%

%% Fit intensity trace

[p,ci] = bleachFit(J,'showfig') 
k_p_fit = p(2); % fitted value for bleaching rate

%% Run kICS

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

%% ICS for guesses

if use_ics_guess
    ics_run = ICSCompiler(J,xi_lags,eta_lags,'subsets',n_ics_subs,...
        'showfigs');
    
    % compute tau = 0 lag in kICS
    n_theta_arr = n_theta*ones(1,length(kSqVector));
    
    % kICS autocorrelation unnormalized
    r_k_tau_0 = kICS3(J-repmat(mean(J,3),[1,1,size(J,3)]),'normByLag','none'); 
    r_k_circ_0 = ellipticInterp(r_k_tau_0(:,:,1),kSqVector,n_theta_arr);
    %
    
    % get param guesses for kICS fitting
    [params_guess,lb,ub] = getkICSGuess(J,ics_run,p,kSqVector,r_k_circ_0,...
        ksq_noise_lb,'showfig');
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

% close all

kICSCorrSubset = mean(abs(r_k_abs),3); % corresponding subset of kICS ACF
% kICSStdDev = std(abs(r_k_abs),0,3)/sqrt(nReps); % standard deviation of ACF (relevant if nReps > 1)

if do_fit
    err = @(params) timeIntkICSFit(params,kSqVectorSubset,tauVector,...
        k_p_fit,T,'normByLag',normByLag,'err',kICSCorrSubset);
	tic
	%
	opts = optimoptions(@fmincon,'Algorithm','interior-point');
	problem = createOptimProblem('fmincon','objective',...
		err,'x0',params_guess,'lb',lb,'ub',ub,'options',opts);
	gs = GlobalSearch; % global search object
	[opt_params,err_min] = run(gs,problem)
	%
	fitTime = toc; disp(['fitTime = ' num2str(fitTime)]); % time to fit
end

%% plotting

ksq2plot = linspace(kSqVectorSubset(1),kSqVectorSubset(end),nPtsFitPlot); % |k|^2 for plotting best fit function result/theory curves

figure()
hold on

color = lines(length(plotTauLags)); 
plotLegend = cell(1,length(plotTauLags));
h_sim_data = zeros(1,length(plotTauLags));
% plot fit and/or theoy curves
for tauInd = 1:length(plotTauLags) % loop and plot over fixed time lag
    if do_fit
        plot(ksq2plot,timeIntkICSFit(opt_params,ksq2plot,plotTauLags(tauInd),k_p_fit,T,'normByLag',normByLag),...
            'Color',color(tauInd,:)) % plot best fit kICS ACF
    end
    if do_theory
        plot(ksq2plot,timeIntkICSFit(sim_info.true_params,ksq2plot,plotTauLags(tauInd),k_p,T,'normByLag',normByLag),...
            '--','Color',color(tauInd,:)); % plot theoretical kICS ACF
    end
end
% plot simulation data
for tauInd = 1:length(plotTauLags)
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
