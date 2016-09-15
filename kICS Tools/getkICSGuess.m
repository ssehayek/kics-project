% getkICSGuess(...) is meant for guessing bounds and parameters prior to
% fitting.
%
function [params_guess,lb,ub] = getkICSGuess(J,ics_run,bleach_params,ksq,r_k_circ_0,ksq_noise_lb,varargin)

figpath = '';

for ii = 1:length(varargin)
    if any(strcmpi(varargin{ii},{'saveFigs'}))
        if ischar(varargin{ii+1})
            figpath = varargin{ii+1};
        else
            warning(['Unknown option for ''',varargin{ii},...
                ''', using default options.'])
        end
    elseif any(strcmpi(varargin{ii},'showFig'))
        show_fig = 1;
    end
end

size_y = size(J,1);
size_x = size(J,2);

sample_area = size_x*size_y;

if ics_run.n_ics_subs == 1
    warning(['Bounds may be inaccurate since ''ICSCompiler.m'' ran',...
        ' without considering multiple subsets. Note the PSF will be set',...
        ' to a single value.'])
end

%% guesses

params_guess = [10*rand rand 2*rand rand 10*rand 1000*rand];

% guess for PSF
params_guess(5) = mean(ics_run.opt_params(:,2));
%

% guess for noise term
density_guess = (mean(ics_run.opt_params(:,1))*pi*mean(ics_run.opt_params(:,2))^2)^(-1);
N_app_guess = density_guess*sample_area;

iksq_noise_lb = find(ksq >= ksq_noise_lb,1,'first');
noise_var_guess = mean(abs(r_k_circ_0(iksq_noise_lb:end)))/sample_area;

eta_p_guess = noise_var_guess/bleach_params(1)^2*N_app_guess/sample_area;
params_guess(6) = eta_p_guess;
%

%% bounds

lb = eps*ones(1,6);
ub = [10,1,2,1,10,Inf];

% bounds for PSF
[lb(5),ub(5)] = deal(min(ics_run.opt_params(:,2)),max(ics_run.opt_params(:,2)));
%

% bounds for noise term
density_lb = (max(ics_run.ub(:,1))*pi*mean(ics_run.opt_params(:,2))^2)^(-1);
N_app_lb = density_lb*sample_area;
eta_p_lb = noise_var_guess/bleach_params(1)^2*N_app_lb/sample_area;
lb(6) = eta_p_lb;

density_ub = (min(ics_run.opt_params(:,1))*pi*mean(ics_run.opt_params(:,2))^2)^(-1);
N_app_ub = density_ub*sample_area;
eta_p_ub = noise_var_guess/bleach_params(1)^2*N_app_ub/sample_area;
ub(6) = eta_p_ub;
%

%% plot and save figures

if ~isempty(figpath) || show_fig
    if ~isempty(figpath) && ~exist(figpath,'dir')
        mkdir(figpath)
    end
    
    figure()
    hold on
    ax = gca;
    
    plot(ksq,abs(r_k_circ_0),'.','markersize',16)
    
    y_min = ax.YLim(1);
    y_max = ax.YLim(2);
    plot([ksq_noise_lb,ksq_noise_lb],[y_min,y_max],'--','linewidth',2)
    
    xlabel('$|\mathbf{k}|^2$ (pixels$^{-2}$)','interpreter','latex','fontsize',14)
    ylabel('$R(|\mathbf{k}|^2,\tau=0)$','interpreter','latex','fontsize',14)
    legend({'data','threshold'},'fontsize',12,'interpreter','latex')
    
    if ~isempty(figpath)
        filename = [figpath,filesep,'noise_thresh','.fig'];
        saveas(gcf,filename);
        filename = [figpath,filesep,'noise_thresh','.pdf'];
        saveas(gcf,filename);
        
        close(gcf)
    end
end
