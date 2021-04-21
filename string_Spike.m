function [strSpikes] = string_Spike(Spikes)
for i = 1:length(Spikes(1,:))
    str = string(Spikes(:,i));
    strSpikes(i,1) = strjoin(str);
end
end

