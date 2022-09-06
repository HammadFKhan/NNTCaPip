function Total_shuffled = allShuffle(Spikes,surrogate)
if nargin<2 || strcmp(surrogate,''), surrogate = 5; end
cell_count = size(Spikes,1);
num_images = size(Spikes,2);
Total_shuffled = zeros(cell_count,num_images);
for x = 1:surrogate
    p = randperm(num_images);
    q = randperm(cell_count)';
    for i = 1:cell_count
        for j = 1:num_images
            Total_shuffled(i,j) = Spikes(q(i),p(j));
        end
    end
end
end

