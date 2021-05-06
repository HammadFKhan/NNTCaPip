function [Image_Stack,num_images,Width,Height] = Image_Reader(filename)

%Obtain values for the image stack that will be used throughout the code
%Save the values in a mat file 'Parameters' for easy loading

ImageNorm = imread(filename);
tiff_info = imfinfo(filename);
num_images = numel(tiff_info);


Width = tiff_info.Width;
Height = tiff_info.Height;


%%Load the filestack and store it in a matrix format


Image_Stack = zeros(Height,Width,num_images);


H = waitbar(0,'Loading Images');
for k = 1:num_images
    waitbar(k/num_images)
    Image = imread(filename,'Index', k);
    IMG = Image(:,:,1);
    Image_Stack(:,:,k) = IMG;
    
end
delete(H)

