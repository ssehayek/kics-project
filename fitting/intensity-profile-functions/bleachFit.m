function [p,ci,varargout] = bleachFit(J,varargin)

alpha = 0.05; % default (1-alpha)*100%=95% confidence intervals
show_fig = 0;
for i = 1:length(varargin)
    if strcmpi(varargin{i},'alpha')
        alpha = varargin{i+1};
    elseif any(strcmpi(varargin{i},'showFig'))
        show_fig = 1;
    end
end

t = 1:size(J,3); % all frame indices
I_t = squeeze(mean(mean(J,1),2)); % mean intensity trace

% bleaching profile for single decay rate model
% p(1): amplitude
% p(2): bleach rate
bleach_profile = @(p,t) p(1)*exp(-p(2)*(t+1)); 

lb = [0,0]; % lower-bound of p
ub = [Inf,1]; % upper-bound
x0 = rand()*[1000,1]; % initial guess

% p = lsqcurvefit(bleach_profile,x0,t',I_t,lb,ub); % LSF
[p,~,resid,~,~,~,jacobian] = lsqcurvefit(bleach_profile,x0,t',I_t,lb,ub); % LSF

% confidence intervals on "p"
% note that "nlparci.m" is supposed to be used with "nlinfit.m" (compatible
% with "lsqcurvefit.m"?)
ci = nlparci(p,resid,'jacobian',jacobian,'alpha',alpha);

% plotting
if show_fig || nargout > 2
    figure()
    hold on
    
    plot(t',I_t) % intensity trace plot
    plot(t',bleach_profile(p,t)) % LSF curve
    
    % labeling
    xlabel('$t$ (frames)','interpreter','latex','fontsize',14)
    ylabel('$\overline{i({\bf r},t)}_t$','interpreter','latex','fontsize',14)
    legend({'simulation','theory'},'fontsize',12,'interpreter','latex')
    
    tightfig(gcf) % no white-space (3rd party package; works for release 2015b)
    
    if nargout > 2
        varargout{1} = gcf;
    end
end