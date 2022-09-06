function [NumActiveNodes,NodeList,NumNodes,NumEdges,SpatialCentroid,SpatialCentroidVariance,ActivityCentroid,ActivityCentroidVariance, ActivityCoords]...
    = Network_Analysis(ROIcentroids,Connected_ROI)
NumActiveNodes=0;
NodeList = 0;
NumNodes=size(ROIcentroids,1);
SpatialCentroid = 'Null';
SpatialCentroidVariance = 'Null';
ActivityCentroid = 'Null';
ActivityCentroidVariance='Null';
SpatialCentroid=mean(ROIcentroids,1);
SpatialCentroidVariance=std(ROIcentroids,0,1);
if ~isempty(Connected_ROI) % checks for valid connected ROI argument
    NumEdges = length(Connected_ROI);
    Cell_1=Connected_ROI(1,1);
    Cell_2=Connected_ROI(1,2);
    NodeList=[Cell_1;Cell_2];
    %Count all Nodes with Edges
    if size(Connected_ROI,1) > 1
        NodeList = vertcat(Connected_ROI(:,1),Connected_ROI(:,2));
        NodeList = unique(NodeList);
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
else
    NumActiveNodes=0;
    NumNodes=0;
    NumEdges=0;
    SpatialCentroid = [];
    SpatialCentroidVariance = [];
    ActivityCentroid = [];
    ActivityCentroidVariance=[];
    ActivityCoords = [];
end
%% Plot structure
% NetworkAnalysis.NumActiveNodes = NumActiveNodes;
% NetworkAnalysis.NodeList = NodeList;
% NetworkAnalysis.NumNodes = NumNodes;
% NetworkAnalysis.NumEdges = NumEdges;
% NetworkAnalysis.SpatialCentroid = SpatialCentroid;
% NetworkAnalysis.SpatialCentroidVariance = SpatialCentroidVariance;
% NetworkAnalysis.ActivityCentroid = ActivityCentroid;
% NetworkAnalysis.ActivityCentroidVariance = ActivityCentroidVariance;



