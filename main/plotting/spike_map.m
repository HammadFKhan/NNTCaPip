function spike_map(DeltaFoverF,time)
% if nargin < 2; time = [];end
% ax = time;
% [grad,~]=colorGradient([0 0 0],[1 1 1],128);
time = 1:(size(DeltaFoverF,2)/30.048);
colormap(jet);
normDelta = zeros(size(DeltaFoverF,1),size(DeltaFoverF,2));
for i = 1:size(DeltaFoverF,1)
    normDelta(i,:) = (DeltaFoverF(i,:)-min(DeltaFoverF(i,:)))/(max(DeltaFoverF(i,:))-min(DeltaFoverF(i,:)));
end
imagesc(normDelta,'Xdata',time);
% set(get(h,'title'),'string','\Delta F/F (%)');
axis on;axis tight;box off;
xlabel('Time (s)'); 
ylabel('Neuron');
end
