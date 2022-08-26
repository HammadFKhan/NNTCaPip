function [DeltaFoverF,dDeltaFoverF,ROI,ROIcentroid,Noise_Power,A] =...
    removeROI(DeltaFoverF,dDeltaFoverF,ROI,ROIcentroid,Noise_Power,A,deleteID)
%% manual curation
% Delete cells based on ID number spefied from contour plots
DeltaFoverF(deleteID,:) = [];
dDeltaFoverF(deleteID,:) = [];
Noise_Power(deleteID,:) = [];
ROI(deleteID) = [];
ROIcentroid(deleteID) = [];
A(:,deleteID) = [];
disp(['ROIs successfully removed!'])
