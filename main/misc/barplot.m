function barplot(data)
%plots std of boxplot from single dimensional data
x = nanmean(data);
bar(nanmean(data)),hold on
y = 1:length(ones(size(x)));
err = (nanstd(data));
buff = [];
scatterOn = 1;
if scatterOn
    for i = 1:size(data,2)
        t = data(:,i);
        buff = [buff;t(t>0)];
        labels = repmat({num2str(i)},length(buff),1);
        scatter(i*ones(length(t(t>0)),1),t(t>0),'filled','jitter','on','jitterAmount',0.1)
    end
end
errorbar(y,x,err)

