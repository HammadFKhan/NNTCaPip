% [filename,path] = uigetfile({'*.tiff';'*.tif'});
% clearvars -except gL23 gL5 week shamL23 shamL5  gL23norm gL5norm
pathname = uigetdir(pwd,'Input Directory');
pathname = fullfile(pathname);
directory = dir(fullfile(pathname,'*.tif'));
L = length(directory);
flag = 0; % Analyze Pser channel
for i = 1:L
    filename = directory(i).name;
    t = Tiff(fullfile(pathname,filename),'r');
    imageData = read(t);
    if flag
        greenChannel = double(imageData(:,:,2)-mean(mean(imageData(:,:,2)))); %normalized Mean subtracted image of green channel
        greenChannelNorm = (greenChannel-min(greenChannel,[],'all'))/(max(greenChannel,[],'all')-min(greenChannel,[],'all'));
        greenChannel1 = mean(greenChannel,1);
        greenChannel2 = mean(greenChannel,2);
        green = mean([greenChannel1 greenChannel2']);
        greenTotal(i) = green;
    else
        redChannel = double(imageData-mean(mean(imageData))); %normalized Mean subtracted image of red channel
        redChannel = (redChannel-min(redChannel,[],'all'))/(max(redChannel,[],'all')-min(redChannel,[],'all'));
        
        redChannel1 = mean(redChannel,1);
        redChannel2 = mean(redChannel,2);
        red = mean([redChannel1 redChannel2']);
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
    end
end
%% 
week = 5;
gL23{week} = redTotal;

gL5{week} = redTotal;
%% 
for i = 4
    figure,barplot(1.5*gL4{i}),ylim([0 30]),box off,set(gca,'TickDir','out'),title(['L4' num2str(i)]);
end
%% erase NAN
for i = 1:5
    buff = gL5{i};
    buff(isnan(buff)) = [];
    gL5new{i} = buff;
    buff2 = gL23{i};
    buff2(isnan(buff2)) = [];
    gL23new{i} = buff2;
end
%% line plot with error
x = 1:5;
gL23new{2} = gL23new{2}/2; gL23new{3} = gL23new{3}/2; 
y = cellfun(@mean,gL23new)*200;
err = (cellfun(@std,gL23new))*100;
figure,errorbar(x,y,err);ylim([0 30]),xlim([0 6]),hold on

x = 1:5;
y = cellfun(@mean,gL5new)*300;
err = (cellfun(@std,gL5new))*100;
errorbar(x,y,err);ylim([0 40]),xlim([0 6])

