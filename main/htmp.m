function h = htmp(corr,ticksize)
if nargin<2 || isempty(ticksize); ticksize = 20; end

h = heatmap(corr);
% map = [0.02 0.631 0.631;
%     1 1 1];
% h.Colormap(map);
h.Colormap = jet;
h.GridVisible = 'off';
h.CellLabelColor = 'none';
h.FontSize = 10;
XLabels = 1:length(h.XDisplayLabels);
CustomXLabels = string(XLabels);
CustomXLabels(mod(XLabels,ticksize) ~= 0) = " ";
CustomYLabels = CustomXLabels;
h.XDisplayLabels = CustomXLabels;
h.YDisplayLabels = CustomYLabels;
h.YDisplayData = flipud(h.YDisplayData);  % equivalent to 'YDir', 'Reverse'
end