function [] = fixAxesScale(ax,varargin)

method = 'none';
for i = 1:2:length(varargin)
    if strcmpi(varargin{i},'method')
        % method used for generating new ticks
        %
        % 'auto': use MATLAB 'auto' option for new ticks and labels
        % 'simple': new ticks from removing every second tick
        % 'none'
        if any(strcmpi(varargin{i+1},{'auto','simple','none'}))
            method = lower(varargin{i+1});
        else
            warning(['Unknown option for ''',varargin{i},...
                ''', using default options.'])
        end
    elseif strcmpi(varargin{i},'axes')
        % choose which plot axes to get new ticks for
        if any(strcmpi(varargin{i+1},{'x','y','both'}))
            ax2change = varargin{i+1};
        end
    end
end

x_tick = ax.XTick;
if iscell(ax.XTickLabel)
    for i = 1:length(x_tick)
        x_tick_label(i) = str2num(ax.XTickLabel{i});
    end
else
    x_tick_label = str2num(ax.XTickLabel);
end

y_tick = ax.YTick;
if iscell(ax.YTickLabel)
    for i = 1:length(y_tick)
        y_tick_label(i) = str2num(ax.YTickLabel{i});
    end
else
    y_tick_label = str2num(ax.YTickLabel);
end

switch method
    case 'auto'
        if any(strcmpi(ax2change,{'x','both'}))
            m_x = (x_tick_label(2)-x_tick_label(1))/(x_tick(2)-x_tick(1));
            b_x = x_tick_label(1)-m_x*x_tick(1);
            
            ax.XTickMode = 'auto';
            ax.XTickLabelMode = 'auto';
            
            % dx_new = ax.XTick(2)-ax.XTick(1);
            % ax.XTick = x_tick(1):dx_new:x_tick(end);
            
            x_tick_new = ax.XTick;
            x_tick_label_new = m_x*x_tick_new + b_x;
            
            ax.XTick = x_tick_new;
            ax.XTickLabel = x_tick_label_new;
        end
        if any(strcmpi(ax2change,{'y','both'}))
            m_y = (y_tick_label(2)-y_tick_label(1))/(y_tick(2)-y_tick(1));
            b_y = y_tick_label(1)-m_y*y_tick(1);
            
            ax.YTickMode = 'auto';
            ax.YTickLabelMode = 'auto';
            
            y_tick_new = ax.YTick;
            y_tick_label_new = m_y*y_tick_new + b_y;
            
            ax.YTick = y_tick_new;
            ax.YTickLabel = y_tick_label_new;
        end
    case 'simple'
        if any(strcmpi(ax2change,{'x','both'}))
            ax.XTick = x_tick(1:2:end);
            ax.XTickLabel = x_tick_label(1:2:end);
        end
        if any(strcmpi(ax2change,{'y','both'}))
            % y-axis
            ax.YTick = y_tick(1:2:end);
            ax.YTickLabel = y_tick_label(1:2:end);
        end
end

% fixes conversion from char array to cell array
if ~iscell(ax.XTickLabel)
    ax.XTickLabel = num2cell(str2num(ax.XTickLabel));
end

if ~iscell(ax.YTickLabel)
    ax.YTickLabel = num2cell(str2num(ax.YTickLabel));
end
