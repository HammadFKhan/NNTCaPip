function Cell_Map_Dice(AverageImage,Connected_ROI,ROIbases,NodeSize,EdgeSize)

if nargin < 6; Edgesize = 5;end
if nargin < 5; NodeSize = 5;end


% Centroid = regionprops(AverageImage,'centroid');
% centroids = cat(1, Centroid.Centroid);

centroids = ROIbases;

[H,W] = size(AverageImage);
% background = ones(H,W);
xlim([0 H])
ylim([0 W])
%%
% imshow(mat2gray(background))
% axis on
% hold on;

[lengthConnected,~] = size(Connected_ROI);
hold on;
radius = NodeSize*ones(size(centroids,1),1);
centerX = centroids(:,1);
centerY = centroids(:,2);
% viscircles([centerY, centerX], radius,'Color','k','LineWidth',1);
% clear centerX centerY 
% count = 1;
% for i = 1:length(centroids)
%     centerX(count,1) = centroids(i,1);
%     centerY(count,1) = centroids(i,2);
%     count = count+1;
% end

% for i = 1:length(centroids)
% %     plot(centroids(i,1),centroids(i,2), 'r+', 'MarkerSize', 10, 'LineWidth', 2);
% 
%     radius = NodeSize; 
%     centerX = centroids(i,1);
%     centerY = centroids(i,2);
% %     circle = viscircles([centerX, centerY], radius);
%     for x = 1:H
%         for y = 1:W
%             dX = abs(x-centerX);
%             dY = abs(y-centerY);
%             pixelrad = (dX)^2 + (dY)^2;
%             if pixelrad <= radius^2
%                 background(x,y) = 0;
%             end
%         end
%     end
% end
color = [0 0 0];
sz = pi*(NodeSize)^2;
% scatter(centerY, centerX,sz,'MarkerFaceColor',color,'MarkerFaceAlpha',1,'MarkerEdgeAlpha',0);


hold off;

% imshow(background)
% map = [0.02 0.631 0.8
%     1 1 1];
% colormap(map);
% hold on
%%
if Connected_ROI(1,:) ~= [0,0,0]
    for j = 1:lengthConnected
        Cell1 = ROIbases(Connected_ROI(j,1),:);
        Cell2 = ROIbases(Connected_ROI(j,2),:);
        corr = Connected_ROI(j,3);
        if corr > 0.7
            LineWidth = EdgeSize;
            Color = [0.02 0.631 0.8];
            x1 = Cell1(1,1);
            x2 = Cell2(1,1);
            y1 = Cell1(1,2);
            y2 = Cell2(1,2);
            %         Color = [0.02 0.631 0.631];
            %         Color = 'k';
            x = line([y1,y2],[x1,x2],'LineWidth',LineWidth,'Color',Color); hold on
            x.Color(4) = 0.6;
        elseif corr >= 0.15 && corr <= 0.7
            LineWidth = EdgeSize;
            Color = [.83 0.45 0.05];
            x1 = Cell1(1,1);
            x2 = Cell2(1,1);
            y1 = Cell1(1,2);
            y2 = Cell2(1,2);
            %         Color = [0.02 0.631 0.631];
            %         Color = 'k';
            x = line([y1,y2],[x1,x2],'LineWidth',LineWidth,'Color',Color); hold on
%             
        end
        %         x1 = Cell1(1,1);
        %         x2 = Cell2(1,1);
        %         y1 = Cell1(1,2);
        %         y2 = Cell2(1,2);
        % %         Color = [0.02 0.631 0.631];
        % %         Color = 'k';
        %         x = line([y1,y2],[x1,x2],'LineWidth',LineWidth,'Color',Color);
        %         x.Color(4) = 0.5;
    end
end
scatter(centerY, centerX,sz,'MarkerFaceColor',color,'MarkerFaceAlpha',1,'MarkerEdgeAlpha',0);