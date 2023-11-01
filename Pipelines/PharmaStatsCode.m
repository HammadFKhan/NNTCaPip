%% Line graph
%Using z score
data = SomOptoM1.Amp(32:64,:);
data = data-mean(data(:,1));
% norm = 2*(data-min(data,[],'all'))/(max(data,[],'all')-min(data,[],'all'))-1;
figure,boxplot(data)
%%
%Using raw power values
data = APVM3.Amp(:,:);
norm = ((data-min(data,[],'all'))/(max(data,[],'all')-min(data,[],'all')));
data = norm-mean(norm(:,1));
figure,boxplot(data),ylim([-1 1]),set(gca,'TickDir','out'),box off

xs = ones(size(data,1),1);
x = [xs,2*xs,3*xs];
figure,plot(x,data,'k.'),hold on,ylim([-1 1]),set(gca,'TickDir','out'),box off
for i = 1:5
 line(x(i,:),data(i,:))
end
%% Beta Event Rate
for i = 1:4
    figure,
    boxplot((SomOptoM1.ER{i})),ylim([0 10]),set(gca,'TickDir','out'),box off
end
%%
%Using raw power values
for i = 1:4
    data(i,:) = mean(cell2mat(compiledNeurons(i).behavioralXcorr),1);
end

figure,boxplot(data),ylim([0 .1]),set(gca,'TickDir','out'),box off

xs = ones(size(data,1),1);
x = [xs,2*xs];
figure,plot(x,data,'k.'),hold on,xlim([0 3]),ylim([0 .1]),set(gca,'TickDir','out'),box off
for i = 1:40
    line(x(i,:),data(i,:))
end