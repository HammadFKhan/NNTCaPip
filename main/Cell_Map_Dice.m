function Cell_Map_Dice(AverageImage,Connected_ROI,ROIbases,NodeSize,EdgeSize)

if nargin < 5; Edgesize = 3;end
if nargin < 4; NodeSize = 6;end


% Centroid = regionprops(AverageImage,'centroid');
% centroids = cat(1, Centroid.Centroid);

centroids = ROIbases;

[H,W] = size(AverageImage);
background = ones(H,W);

%%
% imshow(mat2gray(background))
% axis on
% hold on;

[lengthConnected,~] = size(Connected_ROI);

for i = 1:length(centroids)
%     plot(centroids(i,1),centroids(i,2), 'r+', 'MarkerSize', 10, 'LineWidth', 2);

    radius = NodeSize; 
    centerX = centroids(i,1);
    centerY = centroids(i,2);
%     circle = viscircles([centerX, centerY], radius);
    for x = 1:H
        for y = 1:W
            dX = abs(x-centerX);
            dY = abs(y-centerY);
            pixelrad = (dX)^2 + (dY)^2;
            if pixelrad <= radius^2
                background(x,y) = 0;
            end
        end
    end
end

imshow(background)
map = [0.02 0.631 0.631
    1 1 1];
colormap(map);
hold on
%%


if Connected_ROI(1,:) ~= [0,0,0]
    for j = 1:lengthConnected
        Cell1 = ROIbases(Connected_ROI(j,1),:);
        Cell2 = ROIbases(Connected_ROI(j,2),:);
        corr = Connected_ROI(j,3);
        if corr > 0.5
            LineWidth = 2*EdgeSize;
            Color = [0.05 0.2 0.76];
        elseif corr >= 0.1 && corr <= 0.5
            LineWidth = EdgeSize;
            Color = [0.02 0.631 0.8];
        end
        x1 = Cell1(1,1);
        x2 = Cell2(1,1);
        y1 = Cell1(1,2);
        y2 = Cell2(1,2);
%         Color = [0.02 0.631 0.631];
%         Color = 'k';
        x = line([y1,y2],[x1,x2],'LineWidth',LineWidth,'Color',Color);
        x.Color(4) = 0.4;
    end
end