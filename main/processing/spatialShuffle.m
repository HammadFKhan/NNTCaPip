function Event_shuffled = spatialShuffle(Spikes,surrogate)
cell_count = size(Spikes,1);
num_images = size(Spikes,2);
for i = 1:surrogate
    p = randperm(cell_count)';
    Event_shuffled = zeros(cell_count,num_images);
    for i = 1:num_images
        for j = 1:cell_count
            Event_shuffled(j,i) = Spikes(p(j),i);
        end
    end
end
end
