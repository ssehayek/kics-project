function [] = progBar(n,N,varargin)

% show progress for after every multiple of prog_mult (default: 5%)
prog_mult = floor(N*0.05);
for ii = 1:length(varargin)
    if any(strcmpi(varargin{ii},{'progMult'}))
        if isnumeric(varargin{ii+1}) && varargin{ii+1} <= 1
            prog_mult = varargin{ii+1};
        else
            % use default value for prog_mult
        end
    elseif any(strcmpi(varargin{ii},{'time','stopWatch','keepTime'}))
        if isnumeric(varargin{ii+1})
            time = varargin{ii+1};
        end
    end
end

% if N is small i.e. < 1/floor(prog_mult), show progress at every iteration
prog_mult = max(floor(prog_mult),1);
% only show progress if n is multiple of prog_mult
is_mult = mod(n,prog_mult);
if is_mult == 0
    keep_time = exist('time','var');
    switch keep_time
        case 1
            disp(['progress: ',num2str(n/N*100),'% after ',num2str(time),' s.'])
        otherwise
            disp(['progress: ',num2str(n/N*100),'.'])
    end
end