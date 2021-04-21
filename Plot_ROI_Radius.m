function [hObject,eventdata,handles]=Plot_ROI_Radius(hObject,eventdata,handles)

axes(handles.axes1)
imshow(mat2gray(handles.AverageImage))
hold on
for k=1:length(handles.ROIedges)
    boundary = handles.ROIedges{k};
    plot(boundary(:,2), boundary(:,1),...
    'Color','g','LineWidth',3);
end
axes(handles.axes1)
for k = 1:length(handles.ROIedges)
    boundary = handles.ROIedges{k};
    rndRow = ceil(length(boundary)/(mod(rand*k,7)+1));
    col = boundary(rndRow,2); row = boundary(rndRow,1);
    h = text(col+3, row-3, num2str(handles.AverageImageIndexed(row,col)));
   set(h,'Color','g','FontSize',12,'FontWeight','bold');
end
hold on
[H,W]=size(handles.AverageImage);
circle(H/2,W/2,handles.tolerance)