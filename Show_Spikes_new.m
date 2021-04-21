function spike_plotter_new = Show_Spikes(Spikes)
y = length(Spikes(:,1));
x = length(Spikes(1,:));
spikePlot = zeros(y,x);
baseline = 0;
for i = 1:y
    for ii = 1:x
        if Spikes(i,ii) == 1
            spikePlot(i,ii) = Spikes(i,ii)+baseline;
        else
            spikePlot(i,ii) = -10;
        end
    end
    plot(spikePlot(i,:),'.'); hold on;
    baseline = baseline + 1;
end
ylim ([0,y+1]);
end
    
        
