function corr = correlation_jaccard(Spikes)

[H,~] = size(Spikes);
corr = zeros(H);
for i=1:H
    for j=1:H
        val = jaccard(Spikes(i,:),Spikes(j,:));
        if isempty(val) == 1
            val = 0;
        end
        corr(i,j) = val;
    end
end

for i = 1:H
    if mean(Spikes(i,:)) == 0
        corr(i,:) = 0;
        corr(:,i) = 0;
        corr(i,i) = 1;
    end
end

corr(isnan(corr))=0;