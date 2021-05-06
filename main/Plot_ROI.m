function [hObject,eventdata,handles]=Plot_ROI(hObject,eventdata,handles)

axes(handles.axes1)
hold on
for k=1:length(handles.ROIedges)
    boundary = handles.ROIedges{k};
    plot(boundary(:,2), boundary(:,1),...
    'Color',handles.LineColors,'LineWidth',handles.LineWidth);
end
axes(handles.axes1)
for k = 1:length(handles.ROIedges)
    boundary = handles.ROIedges{k};
    rndRow = ceil(length(boundary)/(mod(rand*k,7)+1));
    col = boundary(rndRow,2); row = boundary(rndRow,1);
    h = text(col+3, row-3, num2str(handles.AverageImageIndexed(row,col)));
   set(h,'Color',handles.FontColors,'FontSize',handles.FontSize,'FontWeight','bold');
end