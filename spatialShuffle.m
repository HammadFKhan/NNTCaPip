function Event_shuffled = spatialShuffle(Spikes,num_images,cell_count)
p = randperm(cell_count)';
Event_shuffled = zeros(cell_count,num_images);
for i = 1:num_images
    for j = 1:cell_count
        Event_shuffled(j,i) = Spikes(p(j),i);
    end
end
end
