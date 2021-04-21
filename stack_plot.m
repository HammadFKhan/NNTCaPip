function stack_plot(DeltaFoverF)

x = length(DeltaFoverF(1,:));
y = length(DeltaFoverF(:,1));
baseline = max(DeltaFoverF,[],'all');
for i = 50:60
    gradient = i/y;
    plot(DeltaFoverF(i,:)+(0.7*baseline),'LineWidth',2,'Color',[.30 .835 .384]); hold on;
    baseline = baseline + max(DeltaFoverF,[],'all');
end

end
