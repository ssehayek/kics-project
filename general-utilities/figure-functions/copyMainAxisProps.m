function [] = copyMainAxisProps(ax_1,ax_2)

ax_2.Title = ax_1.Title;

ax_2.XLabel = ax_1.XLabel;
ax_2.YLabel = ax_1.YLabel;

ax_2.XLim = ax_1.XLim;
ax_2.YLim = ax_1.YLim;

ax_2.XTick = ax_1.XTick;
ax_2.YTick = ax_1.YTick;

ax_2.XTickLabel = ax_1.XTickLabel;
ax_2.YTickLabel = ax_1.YTickLabel;

ax_2.CLim = ax_1.CLim;

if ~isempty(ax_1.Legend)
    legend(ax_2,'show');
end