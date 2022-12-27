function [X] = featureProject(featureData,navgtarget,flag)
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
q(:,2) = V(:,size(featureData,2)-1);
q(:,3) = V(:,size(featureData,2)-2);
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
    scatter3(tProjq1,tProjq2,tProjq3,10,[247 148 29]/255,'filled'); hold on; %[43 57 144]/255
    scatter3(uProjq1,uProjq2,uProjq3,10,[0 0 0]/255,'filled'); %[0 148 68]/255
end

ylabel('bk')
xlabel('ak');
X = [featureData*q(:,1),featureData*q(:,2),featureData*q(:,3)];

end
