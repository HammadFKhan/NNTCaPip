%% ROIcentroidFix
    ROIcentroid = zeros(size(ROI,1),2);
    count = 1;
for i = 1:size(ROI,1)
    if ~isempty(ROI{i})
    ROIcentroid(i,:) = mean(cell2mat(ROI{i}),1);
    end
end
