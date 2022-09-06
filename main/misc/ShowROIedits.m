function ShowROIedits(AverageImage,ROI,L)


imshow(mat2gray(AverageImage)); 
hold on;
colors=['b' 'g' 'r' 'c' 'm' 'y'];
for k=1:length(ROI)
    boundary = ROI{k};
    plot(boundary(:,2), boundary(:,1),...
    'Color','b','LineWidth',3);
end
for k = 1:length(ROI)
    boundary = ROI{k};
    rndRow = ceil(length(boundary)/(mod(rand*k,7)+1));
    col = boundary(rndRow,2); row = boundary(rndRow,1);
    h = text(col+3, row-3, num2str(L(row,col)));
   set(h,'Color','r','FontSize',24,'FontWeight','bold');
end