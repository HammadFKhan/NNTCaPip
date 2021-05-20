function stack_plot(DeltaFoverF,spacing)
if nargin < 2 || strcmp(spacing,'')
    spacing = 1;
end
x = length(DeltaFoverF(1,:));
y = length(DeltaFoverF(:,1));
baseline = spacing*max(DeltaFoverF,[],'all');
for i = 1:y
    gradient = i/y;
    plot(DeltaFoverF(i,:)+(spacing*baseline),'LineWidth',2,'Color',[.30 .835 .384]); hold on;
    baseline = baseline + spacing*max(DeltaFoverF,[],'all');
    axis tight, box off,...
            set(gca,'YTick','','YTickLabel','');
end

end
