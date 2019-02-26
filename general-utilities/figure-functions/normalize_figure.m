function [] = normalize_figure(h,varargin)

%% default parameters

% structure to store all properties
s = struct;
%%% figure properties
% position units of figure
s.units = 'inches';
% width of figure (fixes aspect ratio)
s.width = 3.25;
%%% axis properties
% axis font size
s.fontSize = 9;
% axis title font size
s.titleSize = 11;
% axis label font size
s.labelSize = 11;
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

% reset figure position to default
set(h,'Position',get(0,'defaultfigureposition'));

h.Units = s.units;
% move figure to lower-left of screen
h.Position(1:2) = [0,0];
% change width while keeping aspect ratio
h.Position(3:4) = s.width*h.Position(3:4)/h.Position(3);

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
end

tightfig(h)