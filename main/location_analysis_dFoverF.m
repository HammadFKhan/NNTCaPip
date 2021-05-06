function [Intensity] = location_analysis_dFoverF(ROI, Image_Stack,num_images)

Intensity = zeros(length(ROI),num_images);


H = waitbar(0,'Analyzing High value areas');
for i= 1:length(ROI)
    waitbar(i/length(ROI))
%     Region = cell2mat(ROI{i});
    Region = ROI{i,1};
    [RegionH,~] = size(Region);
    for k = 1:num_images    
        IntensityRegion = 0;
        for j = 1:RegionH
            Area = Region{j,1};
            h = Area(1,1);
            w = Area(1,2);
            Intens = Image_Stack(h,w,k);
            IntensityRegion = (IntensityRegion + Intens);
        end
        val = (IntensityRegion)/RegionH;
        Intensity(i,k) = val;
    end
end
delete(H)


