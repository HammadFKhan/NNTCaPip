function plot_Centroids(SpatialCentroid,SpatialCentroidVariance,ActivityCentroid,ActivityCentroidVariance)



for i = 1:4
    subplot(2,2,i)
    aNetwork=xdev(i,1); % Network Full horizontal radius
    bNetwork=ydev(i,1); % Network Full vertical radius
    x0Network=xcoord(i,1); % x0 Network ellipse centre coordinate
    y0Network=ycoord(i,1); % y0 Network ellipse centre coordinate
    aPhysical=xdev(i,2); % Network Full horizontal radius
    bPhysical=ydev(i,2); % Network Full vertical radius
    x0Physical=xcoord(i,2); % x0 Network ellipse centre coordinate
    y0Physical=ycoord(i,2); % y0 Network ellipse centre coordinate
    %Develop Ellipse Range
    angle=0/180*pi;
    r=0:0.1:2*pi+0.1;
    
    %Newtork Full Plotting
    p=[(aNetwork*cos(r))' (bNetwork*sin(r))'];
    alpha=[cos(angle) -sin(angle)
       sin(angle) cos(angle)];
    p1=p*alpha;
    colorNetwork=[98/255,190/255,112/255];
    patch(x0Network+p1(:,1),y0Network+p1(:,2),colorNetwork,'EdgeColor',colorNetwork,'FaceAlpha',0.75);
    hold on;
    %Nework 1.5 Plotting
    p=[(1.5*aNetwork*cos(r))' (1.5*bNetwork*sin(r))'];
    alpha=[cos(angle) -sin(angle)
       sin(angle) cos(angle)];
    p1=p*alpha;
    colorNetwork=[98/255,190/255,112/255];
    patch(x0Network+p1(:,1),y0Network+p1(:,2),colorNetwork,'EdgeColor',colorNetwork,'EdgeAlpha',0.5,'FaceAlpha',0.5);
    hold on;
    %Physical Full Plotting
    p=[(aPhysical*cos(r))' (bPhysical*sin(r))'];
    alpha=[cos(angle) -sin(angle)
       sin(angle) cos(angle)];
    p1=p*alpha;
    colorPhysical=[120/255,37/255,142/255];
    patch(x0Physical+p1(:,1),y0Physical+p1(:,2),colorPhysical,'EdgeColor',colorPhysical,'FaceAlpha',0.75);
    %Physical 1.5 Plotting
    p=[(1.5*aPhysical*cos(r))' (1.5*bPhysical*sin(r))'];
    alpha=[cos(angle) -sin(angle)
       sin(angle) cos(angle)];
    p1=p*alpha;
    patch(x0Physical+p1(:,1),y0Physical+p1(:,2),colorPhysical,'EdgeColor',colorPhysical,'EdgeAlpha',0.5,'FaceAlpha',0.5);
    
    %Plot Characteristic Edits
    ylim([0,2048])
    xlim([0,2048])
end

%%
for i = 1:4
    subplot(2,2,i)
    aNetwork=xdev(i,1); % Network Full horizontal radius
    bNetwork=ydev(i,1); % Network Full vertical radius
    x0Network=xcoord(i,1); % x0 Network ellipse centre coordinate
    y0Network=ycoord(i,1); % y0 Network ellipse centre coordinate
    aPhysical=xdev(i,2); % Network Full horizontal radius
    bPhysical=ydev(i,2); % Network Full vertical radius
    x0Physical=xcoord(i,2); % x0 Network ellipse centre coordinate
    y0Physical=ycoord(i,2); % y0 Network ellipse centre coordinate
    %Develop Ellipse Range
    angle=0/180*pi;
    r=0:0.1:2*pi+0.1;
    
    %Newtork Full Plotting
    p=[(aNetwork*cos(r))' (bNetwork*sin(r))'];
    alpha=[cos(angle) -sin(angle)
       sin(angle) cos(angle)];
    p1=p*alpha;
    colorNetwork=[98/255,190/255,112/255];
    patch(x0Network+p1(:,1),y0Network+p1(:,2),colorNetwork,'EdgeColor',colorNetwork,'FaceAlpha',0.75);
    hold on;
    %Nework 1.5 Plotting
    p=[(1.5*aNetwork*cos(r))' (1.5*bNetwork*sin(r))'];
    alpha=[cos(angle) -sin(angle)
       sin(angle) cos(angle)];
    p1=p*alpha;
    colorNetwork=[98/255,190/255,112/255];
    patch(x0Network+p1(:,1),y0Network+p1(:,2),colorNetwork,'EdgeColor',colorNetwork,'EdgeAlpha',0.5,'FaceAlpha',0.5);
    hold on;
    %Physical Full Plotting
    p=[(aPhysical*cos(r))' (bPhysical*sin(r))'];
    alpha=[cos(angle) -sin(angle)
       sin(angle) cos(angle)];
    p1=p*alpha;
    colorPhysical=[120/255,37/255,142/255];
    patch(x0Physical+p1(:,1),y0Physical+p1(:,2),colorPhysical,'EdgeColor',colorPhysical,'FaceAlpha',0.75);
    for j=1:10
    %Physical 1.5 Plotting
        p=[((1+0.5*j/10)*aPhysical*cos(r))' ((1+0.5*j/10)*bPhysical*sin(r))'];
        alpha=[cos(angle) -sin(angle)
           sin(angle) cos(angle)];
        p1=p*alpha;
        patch(x0Physical+p1(:,1),y0Physical+p1(:,2),colorPhysical,'EdgeColor',colorPhysical,'EdgeAlpha',0.5,'FaceAlpha',0.5);
    end
    %Plot Characteristic Edits
    ylim([0,2048])
    xlim([0,2048])
end