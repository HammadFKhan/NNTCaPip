function EnsembleMap(AverageImage,ROIbases,NodeList,NodeSize,color)

if nargin < 4; NodeSize = 6;end
if nargin < 5; color = [1 0 0];end
% Centroid = regionprops(AverageImage,'centroid');
% centroids = cat(1, Centroid.Centroid);

centroids = ROIbases;

[H,W] = size(AverageImage);
background = ones(H,W);
xlim([0 H])
ylim([0 W])
%%
% imshow(mat2gray(background))
% axis on
%     plot(centroids(i,1),centroids(i,2), 'r+', 'MarkerSize', 10, 'LineWidth', 2);
hold on;
radius = NodeSize*ones(size(centroids,1),1);
centerX = centroids(:,1);
centerY = centroids(:,2);
viscircles([centerX, centerY], radius,'Color','k','LineWidth',1);
clear centerX centerY 
count = 1;
for i = 1:length(centroids)
    if any(NodeList(:) == i)
    centerX(count,1) = centroids(i,1);
    centerY(count,1) = centroids(i,2);
    count = count+1;
    end
end
if  length(NodeList)>1
    sz = pi*(NodeSize)^2;
    scatter(centerX, centerY,sz,'MarkerFaceColor',color,'MarkerFaceAlpha',1,'MarkerEdgeAlpha',0);
end
hold off;
% imshow(background)
% map = [0.02 0.631 0.8
%     1 1 1];
% colormap(map);
end