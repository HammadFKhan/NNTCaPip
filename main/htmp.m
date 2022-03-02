function h = htmp(corr,ticksize)
if nargin<2 || isempty(ticksize); ticksize = 20; end

h = imagesc(corr);
% [grad,~]=colorGradient([1 1 1],[0 0 0],64);
% colormap(grad);
colormap(jet);
set(get(colorbar,'title'),'string','\Delta F/F (%)');
axis on;axis tight;box off; axis xy;
end