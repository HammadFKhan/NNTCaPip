function h = plot_raster(x,y,width,color,linew)
% Simple XY plot with the symbol '|' instead of '.' or 'o'.

if nargin<1; help plot_raster; return; end
if nargin<2 || isempty(y); y=x; x=1:length(y); end
if nargin<3 || isempty(width); width=1; end
if nargin<4 || isempty(color); color='k'; end
if nargin<5 || isempty(linew); linew=1; end

if length(x)==1 && length(y)>1
    x=ones(size(y))*x;
end
if length(y)==1 && length(x)>1
    y=ones(size(x))*y;
end

plot([x(:)'; x(:)'], [y(:)'+width/2; y(:)'-width/2],'-','color',color,'linewidth',linew)
