function h = htmp(corr,ticksize)
if nargin<2 || isempty(ticksize); ticksize = 20; end

h = imagesc(corr);
% map = [0.02 0.631 0.631;
%     1 1 1];
% h.Colormap(map);
% [grad,~]=colorGradient([.1 .1 .1],[1 1 1],128);
% colormap(grad);
colormap(jet);
set(get(colorbar,'title'),'string','\Delta F/F (%)');
axis on;axis tight;box off; axis xy;