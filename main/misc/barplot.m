function barplot(data)
%plots std of boxplot from single dimensional data
x = mean(data);
bar(mean(data)),hold on
y = 1;
err = (std(data))*ones(size(x));
errorbar(y,x,err)