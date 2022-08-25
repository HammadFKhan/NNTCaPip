%%Define Global Variables
global AverageImage
global ROIcentroid
global ROI
global Image_Stack
global num_images
global DeltaFoverF
global dDeltaFoverF
global Noise_Power
global files
global A
global C
global ops
AverageImage = 0;
ROIcentroid = 0;
ROI = 0;
Image_Stack = 0;
num_images = 0;
DeltaFoverF = 0;
dDeltaFoverF = 0;
Noise_Power = 0;
files = 0;
A = [];
C = [];
ops = [];

% memMap_CaImAn()
CaImAnFull();
