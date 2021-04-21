function [AverageImage] = Average_Image(Image_Stack,num_images,Width,Height,MinVal)

Avg_Image = zeros(Height,Width);




%%
H = waitbar(0,'Creating Average Image');
for k = 1:num_images
    waitbar(k/num_images) 
    for j = 1:Width
        for i = 1:Height
            if Avg_Image(i,j) == 0
                if Image_Stack(i,j,k) <= MinVal
                    Avg_Image(i,j) = 0;
                else
                    Avg_Image(i,j) = 1;
                end
            end
        end
    end
end
delete(H)



%%
Image = imbinarize(Avg_Image(:,:));
Image = bwareaopen(Image(:,:),40);
AverageImage = imfill(Image(:,:),'holes');
end



