rng shuffle % reseed the random number generator; otherwise same random 
            % numbers are generated each time MATLAB restarts

%% Input Variable

%%%%%% Internal Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

codePath = ''; % path where all codes are stored (recursive)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Simulation Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (refer to dronpaSim.m for variable definitions)

nReps = 5; % number of simulations

% main params
sz = 64; 
T = 128; 
w0 = 4; 
N_diff = 500; 
D = 1; 
k_on = 0.1; 
k_off = 0.1; 
k_p = 0; 

% aggregate params
prob_agg = 0;
mean_agg_num = 5;
std_agg_dist = 0.1; 

% filament params 
num_filaments = 20; 
prob_place = 0.3; 

% varargin
snr = 8; % signal-to-noise ratio (enter 'none' for no noise)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% kICS Fitting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

normByLag = 0; % actual lag values i.e. tau=0 is 0th lag

% fitting 
tauVector = 1:5; % tau values to fit

% format [D,k_on,k_off,f_d,w0,sigma]
% here "f_d" is the fraction of diffusing particles, and "sigma" is a term
% which is related to the noise

params_guess = 100*rand(1,6); % guess params for fit
lb = eps*ones(1,6); % lower bound for fit params
ub = [100 1 1 100 100 1000]; % upper bound

% plot kICS ACF
plotTauLags = [1:5]; % tau values to plot
nPtsFitPlot = 1000; % number of |k|^2 pts to plot in best fit function
kSqMin = 0.01; % min |k|^2 value to fit/plot
kSqMax = 3; % max |k|^2 value to fit/plot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Simulations

addpath(genpath(codePath)) % add all codes to path variable (recursive)

J = zeros(sz,sz,T,nReps); 
for n = 1:nReps % create simulations
    [J(:,:,:,n)] = dronpaSim(sz,T,w0,N_diff,D,k_on,k_off,k_p,...
        prob_agg,mean_agg_num,std_agg_dist,num_filaments,prob_place,'snr',snr);
end

%% Run kICS

r_k_norm = zeros(sz,sz,T,nReps);

for n = 1:nReps
    r_k_norm(:,:,:,n) = kICS3(J(:,:,:,n)-repmat(mean(J(:,:,:,n),3),[1,1,T]),'normByLag',normByLag); % kICS autocorrelation function (ACF)
    if n == 1
        [r_k_circ,kSqVector] = circular(r_k_norm(:,:,:,n),floor(sz/2)-2); % circular averaging over same |k|^2 in kICS ACF
        r_k_circ = repmat(r_k_circ,[1,1,nReps]);
    else
        [r_k_circ(:,:,n),~] = circular(r_k_norm(:,:,:,n),floor(sz/2)-2); 
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

kSqFit = linspace(kSqVectorSubset(1),kSqVectorSubset(end),nPtsFitPlot); % |k|^2 for plotting best fit function result

figure()
hold on

color = lines(length(plotTauLags)); 
plotLegend = cell(1,length(plotTauLags));
h_sim_data = zeros(1,length(plotTauLags));
if nReps == 1 % loop and plot over fixed simulation number
    for tauInd = 1:length(plotTauLags) % loop and plot over fixed time lag
        plot(kSqFit,kICSNormTauFitFluctNoise(opt_params,kSqFit,plotTauLags(tauInd),'normByLag',normByLag),'Color',color(tauInd,:))
        h_sim_data(tauInd) = plot(kSqVectorSubset,kICSCorrSubset(:,plotTauLags(tauInd)+1),'.','MarkerSize',16,'Color',color(tauInd,:));
        plotLegend{tauInd} = ['$\tau = ' num2str(plotTauLags(tauInd)) '$'];
    end
else
    for tauInd = 1:length(plotTauLags) % 1:length(plotTauLags)
        plot(kSqFit,kICSNormTauFitFluctNoise(opt_params,kSqFit,plotTauLags(tauInd),'normByLag',normByLag),'Color',color(tauInd,:),'linewidth',1.4)
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
