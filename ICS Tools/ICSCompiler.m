% ICSCompiler(...) is meant for fitting ICS autocorrelation functions
% (ACFs), and to return information about the fit(s).
%
% INPUT PARAMS
% J: the movie to be analyzed (1st 2 dims are spatial, 3rd is temporal)
% xi(eta)_lags: spatial lags to fit in x(y)-direction (note x is dim 2, and
%               y is dim 1.
%

function [ics_run] = ICSCompiler(J,xi_lags,eta_lags,varargin)

% boolean to include (0,0) spatial lag when fitting
includeZeroLag = 0;
% number of subsets to split ACF into (fit is run on each subset)
n_ics_subs = 1;
% path for saving best fits of ICS subsets
figpath = '';
%
show_figs = 0;
for ii = 1:length(varargin)
    if any(strcmpi(varargin{ii},{'includeZeroLag','zeroLag'}))
        if varargin{ii+1} == 0 || varargin{ii+1} == 1
            includeZeroLag = varargin{ii+1};
        else
            warning(['Unknown option for ''',varargin{ii},...
                ''', using default options.'])
        end
    elseif any(strcmpi(varargin{ii},{'nSubsets','subsets'}))
        if varargin{ii+1} >= 1
            n_ics_subs = varargin{ii+1};
        else
            warning(['Unknown option for ''',varargin{ii},...
                ''', using default options.'])
        end
    elseif any(strcmpi(varargin{ii},{'saveFigs'}))
        if ischar(varargin{ii+1})
            figpath = varargin{ii+1};
        else
            warning(['Unknown option for ''',varargin{ii},...
                ''', using default options.'])
        end
    elseif any(strcmpi(varargin{ii},{'showFigs'}))
        show_figs = 1;
    end
end

% compute ICS ACF
[corr,xi_grid,eta_grid] = ICS(J,'meantype','temporal');

% pick out a subset of the ACF
[sub_corr,xi_sub_grid,eta_sub_grid] = getICSCorrSub(corr,xi_lags,...
    eta_lags,'includeZeroLag',includeZeroLag);

% total number of frames
T = size(J,3);
% number of frames in each subset (except maybe last one)
sub_ics_frames = ceil(T/n_ics_subs);

% arrays to store ACF subsets
corr_n = cell(1,n_ics_subs);
sub_corr_n = cell(1,n_ics_subs);
sub_corr_n_avg = zeros(size(sub_corr,1),size(sub_corr,2),n_ics_subs);
for ii = 1:n_ics_subs
    if ii ~= n_ics_subs % not last iteration
        sub_rng = 1+(ii-1)*sub_ics_frames:ii*sub_ics_frames;
    else % last iteration
        sub_rng = 1+(n_ics_subs-1)*sub_ics_frames:T;
    end
    corr_n{:,ii} = corr(:,:,sub_rng);
    sub_corr_n{:,ii} = sub_corr(:,:,sub_rng);
    sub_corr_n_avg(:,:,ii) = mean(sub_corr(:,:,sub_rng),3);
end

% array to store best fit params
opt_params = zeros(n_ics_subs,3);
% array to store objective function values
err_min = zeros(1,n_ics_subs);
% arrays to store initial guesses and bounds
guess_params = zeros(n_ics_subs,3);
lb = zeros(n_ics_subs,3);
ub = zeros(n_ics_subs,3);

disp('fitting ICS subsets...')
tic
parfor ii = 1:n_ics_subs
    corr_ii = corr_n{ii};
    sub_corr_ii = sub_corr_n{ii};
    sub_corr_ii_avg = sub_corr_n_avg(:,:,ii);
    % get initial guess
    [guess_params(ii,:),lb(ii,:),ub(ii,:)] = getICSGuess(corr_ii,sub_corr_ii);
    
    err = @(params) ICSFit(params,xi_sub_grid,eta_sub_grid,'err',sub_corr_ii_avg);
    
    opts = optimoptions(@fmincon,'Algorithm','interior-point');
    problem = createOptimProblem('fmincon','objective',...
        err,'x0',guess_params(ii,:),'lb',lb(ii,:),'ub',ub(ii,:),'options',opts);
    gs = GlobalSearch('Display','off'); % global search object
    [opt_params(ii,:),err_min(ii)] = run(gs,problem);
end
delete(gcp)
toc

% save raw ICS autocorr data
ics_run.corr = corr;
ics_run.xi_grid = xi_grid; ics_run.eta_grid = eta_grid;
ics_run.sub_corr = sub_corr;
ics_run.xi_sub_grid = xi_sub_grid; ics_run.eta_sub_grid = eta_sub_grid;
ics_run.sub_corr_t = sub_corr_n_avg;

% save guess params and bounds
ics_run.guess_params = guess_params;
ics_run.lb = lb;
ics_run.ub = ub;

% save fit results
ics_run.opt_params = opt_params;
ics_run.err_min = err_min;

% save miscellaneous
ics_run.n_ics_subs = n_ics_subs;

%% plot and save figures

if ~isempty(figpath) || show_figs
    if ~isempty(figpath) && ~exist(figpath,'dir')
        mkdir(figpath)
    end
    
    [x_mid_sub,y_mid_sub] = getCtrPxl(sub_corr);
    for ii = 1:n_ics_subs
        figure()
        hold on
        
        plot(eta_lags,sub_corr_n_avg(y_mid_sub,x_mid_sub+eta_lags,ii),'.','markersize',16)
        plot(eta_lags,ICSFit(opt_params(ii,:),0,eta_lags),'-','linewidth',2)
        plot(eta_lags,ICSFit(guess_params(ii,:),0,eta_lags),'--','linewidth',2)
        
        xlabel('$\xi$ (pixels)','fontsize',12,'interpreter','latex')
        ylabel('$\phi(\xi)$','fontsize',12,'interpreter','latex')
        legend({'data','fit','guess'},'fontsize',12,'interpreter','latex')
        
        tightfig(gcf)
        
        if ~isempty(figpath)
            filename = [figpath,filesep,'ics_sub_fit_',num2str(ii),'.fig'];
            saveas(gcf,filename);
            filename = [figpath,filesep,'ics_sub_fit_',num2str(ii),'.pdf'];
            saveas(gcf,filename);
            
            close(gcf)
        end
    end
end
