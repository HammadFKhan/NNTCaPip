function [corr,val] = correlation_dice(Spikes)
fprintf('Calculating Dice correlation...');
[H,~] = size(Spikes);
corr = zeros(H);
for i=1:H
    for j=1:H
        val = dice(Spikes(i,:),Spikes(j,:));
        if isempty(val) == 1
            val = 0;
        end
        corr(i,j) = val;
        if i == j
            corr(i,j) = NaN;
        end
    end
end

for i = 1:H
    if mean(Spikes(i,:)) == 0
        corr(i,:) = 0;
        corr(:,i) = 0;
        corr(i,i) = 1;
    end
end

corr(isnan(corr))=1;

fprintf('done\n')
end