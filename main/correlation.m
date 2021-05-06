function [RHO, PVAL] = correlation(Spikes)

[RHO, PVAL] = corr(Spikes');



[H,~] = size(Spikes);

for i = 1:H
    if mean(Spikes(i,:)) == -10
        PVAL(i,:) = 1;
        PVAL(:,i) = 1;
    end
end

PVAL(isnan(PVAL))=1;
