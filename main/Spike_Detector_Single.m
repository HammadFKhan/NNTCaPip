function Spikes = Spike_Detector_Single(DeltaFoverF,stddev,spikemin)


Spikes = zeros(size(DeltaFoverF,1),size(DeltaFoverF,2),size(DeltaFoverF,3));

% dev = stddev*std(std(DeltaFoverF(:,:,:)));
avghx = 0;%mean(mean(DeltaFoverF(:,:,:)));
for j = 1:size(DeltaFoverF,1)
    tempSpikes = Spikes(j,:);
    dev = stddev*mean(std(DeltaFoverF(j,:,:)));
    temp = DeltaFoverF(j,:);
    tempSpikes(temp>= dev + avghx & temp >= spikemin) = 1;
    Spikes(j,:) = tempSpikes;
%     if DeltaFoverF(j,k) >= dev + avghx && DeltaFoverF(j,k) >= spikemin
%         Spikes(j,k) = 1;
%     else
%     end
end
