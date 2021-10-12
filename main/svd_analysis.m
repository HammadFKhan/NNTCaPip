%% SVD test
function svd_analysis
figure,subplot(2,2,1), axis off
h = htmp(Wi);
plotind = 2;
X = Wi;
[U,S,V] = svd(X);
for r = [2 3 5]
Xapprox = U(:,r)*S(r,r)*V(:,r)';
subplot(2,2,plotind), plotind = plotind+1;
htmp(Xapprox);
caxis([0 1]);
end
%% Mean-substracted data
% figure,subplot(2,2,1), axis off
% h = htmp(sim_index);h.ColorbarVisible = 'off';
% plotind = 2;
% avgSim = mean(sim_index,2);
% X = sim_index - avgSim*ones(1,size(sim_index,2));
% [U,S,V] = svd(X);
% for r = [2 3 4]
% Xapprox = U(:,r)*S(r,r)*V(:,r)';
% subplot(2,2,plotind), plotind = plotind+1;
% h = htmp(Xapprox);
% h.ColorbarVisible = 'off';caxis([0 0.15]);saveas(gcf,'myName','png');
% end

figure,subplot(2,2,1), axis off
h = htmp(sim_index);caxis([0 1]);h.ColorbarVisible = 'off';
plotind = 2;
avgSim = mean(sim_index,2);
X = sim_index - avgSim*ones(1,size(sim_index,2));
[U,S,V] = svd(X);
for r = [1 3 5]
Xapprox = U(:,r)*S(r,r)*V(:,r)';
subplot(2,2,plotind), plotind = plotind+1;
h = htmp(Xapprox);caxis([0 1]);
h.ColorbarVisible = 'off';caxis([0 0.15]);h.Colormap = flipud(gray(5));saveas(gcf,'myName','png');
end
%% truncate 40x40

figure,subplot(2,2,1), axis off
sim_index2 = sim_index1(1:40,1:40);
h = htmp(sim_index2);h.ColorbarVisible = 'off';
plotind = 2;
avgSim2 = mean(sim_index2,2);
X = sim_index2 - avgSim2*ones(1,size(sim_index2,2));
[U,S,V] = svd(X);
for r = [1 2 3 4 5 6 7 8 9 10]
Xapprox = U(:,r)*S(r,r)*V(:,r)';
subplot(2,2,plotind), plotind = plotind+1;
h = htmp(Xapprox);
h.ColorbarVisible = 'off';caxis([0 0.15]);h.Colormap = flipud(gray(5));saveas(gcf,'myName','png');
end

%% Singular Values
figure,subplot(1,2,1)
semilogy(diag(S),'k','Linewidth',2),grid on;hold on;
xlabel('r'); ylabel('Singular Value, \Sigma_r')
set(gca,'FontSize',14);
subplot(1,2,2)
plot(cumsum(diag(S))/sum(diag(S)),'k','LineWidth',2),grid on;
xlabel('r'); ylabel('Cumulative Energy');
set(gca,'FontSize',14);

%% Alpha Projection (PCA)
figure,hold on
for i = 1:size(X,1)
    x(i) = V(:,1)'*X(i,:)';
    y(i) = V(:,2)'*X(i,:)';
    z(i) = V(:,3)'*X(i,:)';
end
sz = 10; 
scatter3(x,y,z,sz,'MarkerEdgeColor',[0 .5 .5],...
  'MarkerFaceColor',[0 .7 .7],...
  'LineWidth',1.5),view(-50,50);
xlabel('PC1'),ylabel('PC2'),zlabel('PC3');
grid on,axis([-.5 .5 -.5 .5 -.5 .5]);
%%
figure,

hold on
for j = 1:58
    for i = 1:size(X,1)
        x(i) = V(:,j)'*X(i,:)';
        y(i) = V(:,j+1)'*X(i,:)';
        z(i) = V(:,j+2)'*X(i,:)';
    end
    sz = 5;
    scatter3(x,y,z,sz,'filled'),view(-50,50);
end
xlabel('PC1'),ylabel('PC2'),zlabel('PC3');
axis([-.5 .5 -.5 .5 -.5 .5]);
%%
PCX = x';
PCY = y';
PCZ = z';
%%
figure
for iii = [2,3]
    sz = 48; 
    gradient = 0+iii/3;
    scatter(PCX(:,iii),PCY(:,iii),sz,'MarkerEdgeColor',[1 1 1],...
      'MarkerFaceColor',[gradient 1-gradient .6],'LineWidth',.0001);hold on;
    xlabel('PC1'),ylabel('PC2');
    grid on;
    axis([-1 1 -1 1]);
end
% legend('PFF-Cor','PFF-Str','WT');
set(gca,'FontSize',20)
% set(gca,'YTick',[]);set(gca,'XTick',[]);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);...
    print(gcf,'-painters','-depsc', 'PCA.eps', '-r250'); 



%% Shuffle

[vectorized,sim_index] = cosine_similarity(Spikes,20);
avgSim = mean(sim_index,2);
X = sim_index - avgSim*ones(1,size(sim_index,2));
[U,S,V] = svd(X);
figure,plot(cumsum(diag(S))/sum(diag(S)),'k','LineWidth',2),grid on;hold on;
xlabel('PC #'); ylabel('Cumulative Energy');
set(gca,'FontSize',14);
k = 1;
for j = [10 20 40 100]
    surrogate = j;
    k = k+1;
    gradient = k/6;
    Total_shuffled = allShuffle(Spikes,num_images,cell_count,surrogate);
    [shufvectorized,shufsim_index] = cosine_similarity(Total_shuffled,20);

    avgSim = mean(shufsim_index,2);
    X = shufsim_index - avgSim*ones(1,size(shufsim_index,2));
    [U,S,V] = svd(X);
    plot(cumsum(diag(S))/sum(diag(S)),'Color',[gradient 1-gradient .6],'LineWidth',2),grid on;
    xlabel('r'); ylabel('Cumulative Energy');
    set(gca,'FontSize',14);box off;
end

%% Trial SVD/PCA Analysis
% trial_sim_index = zeros(length(batchData(1).sim_index(:)),length(batchData));
% for jj = 1:length(batchData)
%     trial_sim_index(:,jj) = batchData(jj).sim_index(:);
% end
% figure, axis off
% h = htmp(sim_index);caxis([0 1]);h.ColorbarVisible = 'off';
% plotind = 2;
vectorized = batchData(1).vectorized;
trial_sim_index = vectorized;
avgSim = mean(trial_sim_index,2);
X = trial_sim_index - avgSim*ones(1,size(trial_sim_index,2));
[U,S,V] = svd(X);

% figure,hold on
% for i = 1:size(X,1)
%     x(i) = V(:,1)'*X(i,:)';
%     y(i) = V(:,2)'*X(i,:)';
%     z(i) = V(:,3)'*X(i,:)';
% end
% sz = 10; 
% scatter3(x,y,z,sz,'MarkerEdgeColor',[0 .5 .5],...
%   'MarkerFaceColor',[0 .7 .7],...
%   'LineWidth',1.5),view(-50,50);
% xlabel('PC1'),ylabel('PC2'),zlabel('PC3');
% grid on,axis([-.5 .5 -.5 .5 -.5 .5]);
figure
plot(cumsum(diag(S))/sum(diag(S)),'Color','k','LineWidth',2),grid on;
xlabel('PC #'); ylabel('Cumulative Energy');
set(gca,'FontSize',14);box off;

%% PCA Analysis
[coeff, score,~,~,explained] = pca(batchData(1).sim_index);
figure(),plot(cumsum(explained))
%%
X = sim_index;
XXt = X*X';
XtX = X'*X;
figure,subplot(1,2,1)
imagesc(XtX); axis off;
subplot(1,2,2)
imagesc(XXt);axis off;
end
