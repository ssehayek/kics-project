rng shuffle % reseed the random number generator; otherwise same random 
            % numbers are generated each time MATLAB restarts

%% Input Variable

%%%%%% Internal Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

codePath = 'C:\Users\SimonS\Dropbox (Personal)\Research\PhD\SOFI-Project'; % path where all codes are stored (recursive)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Simulation Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (refer to dronpaSim.m for variable definitions)

nReps = 1; % number of simulations

% main params
sz = 64; 
T = 2048; 
w0 = 4; 
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

%%%%%% kICS Fitting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

do_fit = 0; % run fit algorithm
do_theory = 1; % superimpose theory curves on simulation curves

normByLag = 0; % actual lag values i.e. tau=0 is 0th lag

% fitting 
tauVector = 1:5; % tau values to fit

% format [D,k_on,k_off,f_d,w0,eta']
% here "f_d" is the fraction of diffusing particles, and "sigma" is a term
% which is related to the noise

params_guess = 100*rand(1,6); % guess params for fit
lb = eps*ones(1,6); % lower bound for fit params
ub = [100 1 1 1 100 1000]; % upper bound

% plot kICS ACF
plotTauLags = [1:5]; % tau values to plot
nPtsFitPlot = 1000; % number of |k|^2 pts to plot in best fit function
kSqMin = 0.01; % min |k|^2 value to fit/plot
kSqMax = 3; % max |k|^2 value to fit/plot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Simulations

addpath(genpath(codePath)) % add all codes to path variable (recursive)

% storing true params
true_params = zeros(nReps,6); % array to store true params of simulation
%

J = zeros(sz,sz,T,nReps); 
for n = 1:nReps % create simulations
    [J(:,:,:,n),sim_info] = dronpaSim(sz,T,w0,N_diff,D,k_on,k_off,k_p,...
        prob_agg,mean_agg_num,std_agg_dist,num_filaments,prob_place,'snr',snr);

    true_params(n,:) = sim_info.true_params; % determine true params for each simulation
end
true_params_mean = mean(true_params,1); % mean of true params over all simulations
                                        % note only f_d and sigma are
                                        % variable across simulations

%% Fit intensity trace

t = 1:T;
I_t = squeeze(mean(mean(J,1),2)); % mean intensity trace

bleach_profile = @(x,t) x(1)*exp(-x(2)*t); % bleaching profile for single decay rate model
% x(1): amplitude
% x(2): bleach rate
lb_bleach = [0,0]; % lower-bound of x
ub_bleach = [Inf,1]; % upper-bound
x0 = [1,rand()]; % initial guess

x = lsqcurvefit(bleach_profile,x0,t',I_t,lb_bleach,ub_bleach) % LSF
k_p_fit = x(2); % fit value for k_p

figure()
hold on

plot(t',I_t)
plot(t',bleach_profile(x,t))

%% Run kICS

r_k_norm = zeros(sz,sz,T,nReps);

for n = 1:nReps
    r_k_norm(:,:,:,n) = kICS3(J(:,:,:,n)-repmat(mean(J(:,:,:,n),3),[1,1,T]),'normByLag',normByLag); % kICS autocorrelation function (ACF)
    if n == 1
        [r_k_circ,kSqVector] = circular(r_k_norm(:,:,:,n)); % circular averaging over same |k|^2 in kICS ACF
        r_k_circ = repmat(r_k_circ,[1,1,nReps]);
    else
        [r_k_circ(:,:,n),~] = circular(r_k_norm(:,:,:,n));
    end
    r_k_abs = abs(r_k_circ); % get rid of complex values
end
%% kICS Corr Plot

% close all

kSqMinIndex = find(kSqVector >= kSqMin,1,'first') % lowest index, i, which satisfies kSqVector(i) >= kSqMin
kSqMaxIndex = find(kSqVector <= kSqMax,1,'last') % highest index, j, which satisfies kSqVector(j) <= kSqMax
kSqVectorSubset = kSqVector(kSqMinIndex:kSqMaxIndex); % all values satisfying kSqMin <= kSqVector <= kSqMax
kSqSubsetInd = kSqMinIndex:kSqMaxIndex; % all indices satisfying kSqMin <= kSqVector(i) <= kSqMax

kICSCorrSubset = mean(abs(r_k_abs(kSqSubsetInd,:,:)),3); % corresponding subset of kICS ACF
kICSStdDev = std(abs(r_k_abs(kSqSubsetInd,:,:)),0,3)/sqrt(nReps); % standard deviation of ACF (relevant if nReps > 1)

if do_fit
	err = @(params) kICSNormTauFitFluctNoise(params,kSqVectorSubset,tauVector,'normByLag',normByLag,'err',kICSCorrSubset); % function handle of LSF error
	
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

kSq2Plot = linspace(kSqVectorSubset(1),kSqVectorSubset(end),nPtsFitPlot); % |k|^2 for plotting best fit function result/theory curves

figure()
hold on

color = lines(length(plotTauLags)); 
plotLegend = cell(1,length(plotTauLags));
h_sim_data = zeros(1,length(plotTauLags));
% plot fit and/or theoy curves
for tauInd = 1:length(plotTauLags) % loop and plot over fixed time lag
    if do_fit
        plot(kSq2Plot,kICSNormTauFitFluctNoise(opt_params,kSq2Plot,plotTauLags(tauInd),'normByLag',normByLag),...
            'Color',color(tauInd,:),'linewidth',1.4)
    elseif do_theory
        if k_p == 0
            plot(kSq2Plot,kICSNormTauFitFluctNoise(true_params_mean,kSq2Plot,plotTauLags(tauInd),'normByLag',normByLag),...
                '--','Color',color(tauInd,:),'linewidth',1.4)
        else
            plot(kSq2Plot,kICSNormTauFitFluctNoiseBleach(true_params_mean,kSq2Plot,plotTauLags(tauInd),k_p,T,'normByLag',normByLag),...
                '--','Color',color(tauInd,:),'linewidth',1.4)
        end
    end
end
% plot simulation data
if nReps == 1 
    for tauInd = 1:length(plotTauLags) 
		h_sim_data(tauInd) = plot(kSqVectorSubset,kICSCorrSubset(:,plotTauLags(tauInd)+1),'.','MarkerSize',16,'Color',color(tauInd,:));
        plotLegend{tauInd} = ['$\tau = ' num2str(plotTauLags(tauInd)) '$'];
    end
else
    for tauInd = 1:length(plotTauLags) 
        h_sim_data(tauInd) = errorbar(kSqVectorSubset,kICSCorrSubset(:,plotTauLags(tauInd)+1),...
            kICSStdDev(:,plotTauLags(tauInd)+1),'.','Color',color(tauInd,:),'MarkerSize',12);
        plotLegend{tauInd} = ['$\tau = ' num2str(plotTauLags(tauInd)) '$']; 
    end    
end
% labeling
xlabel('$|\mathbf{k}|^2$ (pixels$^{-2}$)','interpreter','latex','fontsize',14)
ylabel('$\phi(|\mathbf{k}|^2,\tau)$','interpreter','latex','fontsize',14)
legend(h_sim_data,plotLegend,'fontsize',12,'interpreter','latex')
xlim([kSqMin kSqMax])
ylims = get(gca,'ylim');
ylim([0 ylims(2)])

tightfig(gcf)
