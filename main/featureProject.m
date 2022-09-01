function [tProjq1, tProjq2, uProjq1, uProjq2] = featureProject(featureData,navgtarget,flag)
featureData = double(featureData');
covmatrix = (featureData'*featureData);
covmatrix = covmatrix/size(featureData,1);
figure();
imagesc(covmatrix);
colormap(jet);
colorbar;
[V,D] = eig(covmatrix);
figure,
semilogx(flip(diag(D)),'k','Linewidth',2),grid on;hold on;
q(:,1) = V(:,size(featureData,2));
q(:,2) = V(:,size(featureData,2)-2);
q(:,3) = V(:,size(featureData,2)-3);
figure();
plot(q);
ylabel('Voltage (\mu V)')
xlabel('Time');

tProjq1 = featureData(1:navgtarget,:)*q(:,1);
tProjq2 = featureData(1:navgtarget,:)*q(:,2);
tProjq3 = featureData(1:navgtarget,:)*q(:,3);
uProjq1 = featureData(navgtarget+1:end,:)*q(:,1);
uProjq2 = featureData(navgtarget+1:end,:)*q(:,2);
uProjq3 = featureData(navgtarget+1:end,:)*q(:,3);
if flag
    figure()
    scatter(tProjq1,tProjq2,200,'b.'); hold on;
    scatter(uProjq1,uProjq2,200,'r.');
else
    scatter3(tProjq1,tProjq2,tProjq3,80,[43 57 144]/255,'filled'); hold on;
    scatter3(uProjq1,uProjq2,uProjq3,80,[0 148 68]/255,'filled');
end

ylabel('bk')
xlabel('ak');
end
