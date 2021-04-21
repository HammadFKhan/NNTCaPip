function Spikes_shuffled = tempShuffle(Spikes,num_images,cell_count)
p = randperm(num_images);
Spikes_shuffled = zeros(cell_count,num_images);
for i = 1:cell_count
    for j = 1:num_images
        Spikes_shuffled(i,j) = Spikes(i,p(j));
    end
end
end

