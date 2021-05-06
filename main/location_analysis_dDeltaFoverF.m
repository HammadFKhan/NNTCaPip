function [DeltaFoverF] = location_analysis_dDeltaFoverF(ROI, Image_Stack,num_images)

% clear
% clc
% 
% load('Incubator Paper/nNP + Heat/Data/FinalROI Data.mat')
% %%
% %(ROI(:,1), Image_StackH0,num_imagesH0);
% 
% Image_Stack = Image_StackH0;
% num_images = num_imagesH0;
% 
% ROI = ROI(:,1);

%%
AvgIntensityRegion = zeros(length(ROI),num_images);
DeltaFoverF = zeros(length(ROI),num_images-1);

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
        AvgIntensityRegion(i,k) = val;
    end
    for l = 1:num_images-1
          DeltaFoverF(i,l) = (AvgIntensityRegion(i,l+1)-AvgIntensityRegion(i,l))/mean(AvgIntensityRegion(i,:));
%          DeltaFoverF(i,l) = (AvgIntensityRegion(i,l+1)-AvgIntensityRegion(i,l));
    end
end
delete(H)





% deltaI = zeros(length(ROI),num_images-1);
% for ii = 1:length(ROI)
%     for jj = 1:num_images-1
%         deltaI(ii,jj) = AvgIntensityRegion(ii,jj+1) - AvgIntensityRegion(ii,jj);
%     end
% end
