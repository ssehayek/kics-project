function [] = normalize_figure(h,varargin)

%% default parameters

% structure to store all properties
s = struct;
%%% figure properties
% position units of figure
s.units = 'inches';
% boolean for adjusting width
s.changeFigPos = 1;
% width of figure
s.width = 3.25;
%%% axis properties
% axis font size
s.fontSize = 9;
% axis title font size
s.titleSize = 9;
% axis label font size
s.labelSize = 9;
% axis legend font size
s.legSize = 9;
% box
s.box = 'on';

fields = fieldnames(s);
for iVar = 1:length(varargin)
    if any(strcmp(varargin{iVar},fields))
        s.(varargin{iVar}) = varargin{iVar+1};
    end
end

%% 

% move figure to lower-left of screen
if s.changeFigPos
    % reset figure position to default
    set(h,'Position',get(0,'defaultfigureposition'));
    
    h.Units = s.units;
    h.Position(1:2) = [0,0];
    % change width while keeping aspect ratio
    h.Position(3:4) = s.width*h.Position(3:4)/h.Position(3);
end
h.PaperPositionMode = 'auto';
paper_pos = h.PaperPosition;
h.PaperSize = [paper_pos(3) paper_pos(4)];

AxList = findall(h,'type','axes');
for iAx = 1:length(AxList)
    ax = AxList(iAx);
    
    ax.FontSize = s.fontSize;
    ax.Title.FontSize = s.titleSize;
    ax.XLabel.FontSize = s.labelSize;
    ax.YLabel.FontSize = s.labelSize;
    if ~isempty(ax.Legend)
        ax.Legend.FontSize = s.legSize;
    end
    ax.Box = s.box;
    % fixes conversion from char array to cell array
%     fixAxesScale(ax);
end

% tightfig(h)