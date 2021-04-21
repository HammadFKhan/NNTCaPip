function [NumActiveNodes,NumNodes,NumEdges,SpatialCentroid,SpatialCentroidVariance,ActivityCentroid,ActivityCentroidVariance] = Network_Analysis(ROIcentroids,Connected_ROI)
NumActiveNodes=0;
NumNodes=size(ROIcentroids,1);
SpatialCentroid = 'Null';
SpatialCentroidVariance = 'Null';
ActivityCentroid = 'Null';
ActivityCentroidVariance='Null';
SpatialCentroid=mean(ROIcentroids,1);
SpatialCentroidVariance=std(ROIcentroids,0,1);
if Connected_ROI(1,3) ~= 0
    NumEdges = length(Connected_ROI);
    Cell_1=Connected_ROI(1,1);
    Cell_2=Connected_ROI(1,2);
    NodeList=[Cell_1,Cell_2];
    index=3;
    %Count all Nodes with Edges
    if size(Connected_ROI,1) > 1
        for j = 2:size(Connected_ROI,1)
            Cell_1=Connected_ROI(j,1);
            Cell_2=Connected_ROI(j,2);
            if ismember(Cell_1,NodeList)== 0 
                NodeList(index)=Cell_1;
                index=index+1;
            elseif ismember(Cell_2,NodeList)==0
                NodeList(index)=Cell_2;
                index=index+1;
            end
        end
        NumActiveNodes=length(NodeList);
    %Find Centroid of Edges
        ActivityCoords=zeros(size(Connected_ROI,1),2);
         for i=1:size(Connected_ROI,1)
            coord_Cell_1=ROIcentroids(Connected_ROI(i,1),:);
            coord_Cell_2=ROIcentroids(Connected_ROI(i,2),:);
            ActivityCoords(i,:) = (coord_Cell_1+coord_Cell_2)./2;
         end
        ActivityCentroid = mean(ActivityCoords,1);
        ActivityCentroidVariance=std(ActivityCoords,0,1);
    end
else
    NumActiveNodes=0;
    NumNodes=0;
    NumEdges=0;
    SpatialCentroid = 'Null';
    SpatialCentroidVariance = 'Null';
    ActivityCentroid = 'Null';
    ActivityCentroidVariance='Null';
end
