function pearson = correlation_pearson(Spikes)

[H,~] = size(Spikes);
pearson = zeros(H);
for i=1:H
%     for j=1:H
        R = corrcoef(Spikes(i,:),Spikes(i,:));
%         if isempty(R) == 1
%             R = 0;
%         end
%         pearson(i,j) = R;
%     end
end

for i = 1:H
    if mean(Spikes(i,:)) == 0
        pearson(i,:) = 0;
        pearson(:,i) = 0;
        pearson(i,i) = 1;
    end
end

pearson(isnan(pearson))=0;