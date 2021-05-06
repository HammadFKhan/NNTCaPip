function rateImage = firing_rate(Spikes,height,time)

AvgSpikeperFrame = mean(Spikes); %Average network activity per frame
AvgSpikeperROI = mean(Spikes,2); %Average spike rate per ROI

x = length(AvgSpikeperFrame); % Extract x value
y = AvgSpikeperFrame; % Extract rate value

minX = min(x);
maxX = max(x);
minY = min(y);
maxY = max(y);
% height = 30;
rateImage = zeros(height,x);
for i = 1:height
        spikeImage(i,:) = y; 
end
imshow(spikeImage,'Xdata',time);
axis on;axis tight;box off;
h = colorbar; caxis([0 ceil(maxY*10)/10]);
set(get(h,'title'),'string','Firing Rate (FR)');
xlabel('Time (s)');
ax = gca;ax.YAxis.Visible = 'off';
% Colormap is not gray scale.
% Apply some other colormap if you want
colormap(jet(256));
end
