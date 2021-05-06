function ShowROI(AverageImage,ROI,L)


imshow(mat2gray(AverageImage)); 
hold on;
colors=['b' 'g' 'r' 'c' 'm' 'y'];
for k=1:length(ROI)
    boundary = ROI{k}{1,1};
    cidx = mod(k,length(colors))+1;
    plot(boundary(:,2), boundary(:,1),...
    colors(cidx),'LineWidth',2);
    rndRow = ceil(length(boundary)/(mod(rand*k,7)+1));
    col = boundary(rndRow,2); row = boundary(rndRow,1);
    h = text(col+1, row-1, num2str(L(row,col)));
    set(h,'Color',colors(cidx),'FontSize',14,'FontWeight','bold');
end
