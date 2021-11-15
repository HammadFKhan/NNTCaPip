function Spikes_shuffled = tempShuffle(Spikes,surrogate)
cell_count = size(Spikes,1);
num_images = size(Spikes,2);
Spikes_shuffled = Spikes;
for ii = 1:surrogate
    disp(['Surrogate #: ' num2str(ii)])
    for i = 1:cell_count
        p = randperm(num_images);
        Spikes_shuffled(i,:) = Spikes_shuffled(i,p);
    end
end

