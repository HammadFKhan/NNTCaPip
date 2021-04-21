function [coactive_cells,detected_spikes] = coactive_index(Spikes)
for i = 1:length(Spikes(1,:))
    total = length(Spikes(:,i));
    detected_spikes(i) = sum(Spikes(:,i));
    coactive_cells(:,i) = detected_spikes(i)./total;
end
coactive_cells = coactive_cells.*100;
end

    