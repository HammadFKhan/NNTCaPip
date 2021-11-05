function Spikes_shuffled = tempShuffle(Spikes,surrogate)
for ii = 1:surrogate
cell_count = size(Spikes,1);
num_images = size(Spikes,2);
Spikes_shuffled = zeros(cell_count,num_images);
for i = 1:cell_count
    p = randperm(num_images);
    Spikes_shuffled(i,:) = Spikes(i,p);
end
end

