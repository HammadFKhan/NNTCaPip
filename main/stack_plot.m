function stack_plot(DeltaFoverF,spacing,amplitude)
if nargin < 3 || strcmp(amplitude,'')
    amplitude = 1;
end
if nargin < 2 || strcmp(spacing,'')
    spacing = 1;
end
x = length(DeltaFoverF(1,:));
y = length(DeltaFoverF(:,1));
baseline = spacing*max(DeltaFoverF,[],'all');
for i = 1:y
    gradient = i/y;
    plot(amplitude*DeltaFoverF(i,:)+(spacing*baseline),'LineWidth',1,'Color',[.30 .835 .384]); hold on;
    baseline = baseline + spacing*max(DeltaFoverF,[],'all');
    axis tight, box off,...
            set(gca,'YTick','','YTickLabel','');
end

end
