function Cell_Map_Black(AverageImage,Connected_ROI,ROIbases)

Centroid = regionprops(AverageImage,'centroid');
centroids = cat(1, Centroid.Centroid);


[H,W] = size(AverageImage);

background = ones(H,W);


%%
% imshow(mat2gray(background))
% axis on
% hold on;

lengthConnected = length(Connected_ROI);

for i = 1:length(centroids)
%     plot(centroids(i,1),centroids(i,2), 'r+', 'MarkerSize', 10, 'LineWidth', 2);

    radius = 5; 
    centerX = centroids(i,2);
    centerY = centroids(i,1);

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
hold on
%%


if Connected_ROI(1,:) ~= [0,0,1]
    for j = 1:lengthConnected
        Cell1 = ROIbases(Connected_ROI(j,1),:);
        Cell2 = ROIbases(Connected_ROI(j,2),:);
        x1 = Cell1(1,2);
        x2 = Cell2(1,2);
        y1 = Cell1(1,1);
        y2 = Cell2(1,1);
        Color = 'b';
        x = line([y1,y2],[x1,x2],'LineWidth',4,'Color',Color);
        x.Color(4) = 0.3;
    end
end


