% [filename,path] = uigetfile({'*.tiff';'*.tif'});
clearvars -except gL23 gL5 week shamL23 shamL5  gL23norm gL5norm
pathname = uigetdir(pwd,'Input Directory');
pathname = fullfile(pathname);
directory = dir(fullfile(pathname,'*.tif'));
L = length(directory);
for i = 1:L
    filename = directory(i).name;
    t = Tiff(fullfile(pathname,filename),'r');
    imageData = read(t);
    redChannel = double(imageData(:,:,1)-mean(mean(imageData(:,:,1)))); %normalized Mean subtracted image of red channel
%     greenChannel = double(imageData(:,:,2)-mean(mean(imageData(:,:,2)))); %normalized Mean subtracted image of green channel
    redChannel = (redChannel-min(redChannel,[],'all'))/(max(redChannel,[],'all')-min(redChannel,[],'all'));
%     greenChannelNorm = (greenChannel-min(greenChannel,[],'all'))/(max(greenChannel,[],'all')-min(greenChannel,[],'all'));
    
    
    redChannel1 = mean(redChannel,1);
    redChannel2 = mean(redChannel,2);
%     greenChannel1 = mean(greenChannel,1);
%     greenChannel2 = mean(greenChannel,2);
    
    
    red = mean([redChannel1 redChannel2']);
%     green = mean([greenChannel1 greenChannel2']);
    % %%
    % figure,
    % subplot(1,3,1),imagesc(redChannel)
    % subplot(1,3,2),imagesc(greenChannel)
    % subplot(1,3,3),imagesc(imageData)
    % %%
    % figure,hold on
    % area(smoothdata(red,'gaussian',5))
    % area(smoothdata(green,'gaussian',5))
    % legend('Red Channel','Green Channel')
    redTotal(i) = red; 
% greenTotal(i) = green;
end
%% 
week = 3;
gL23{week} = redTotal;

gL5{week} = redTotal;
