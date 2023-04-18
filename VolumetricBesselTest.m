pathname = uigetdir(pwd,'Input Directory');
pathname = fullfile(pathname);
directory = dir(fullfile(pathname,'*.mat'));
L = length(directory);
%%
for i = 1:L
    disp(['Parsing: ' num2str(directory(i).name)])
    S(i).FOV = load(fullfile(directory(i).folder,directory(i).name));
end 
%%
ROIcentroidL1 = [];ROInew = [];
for i = 1:length(dataL1.ROI)
    trans = rand(1)*350;
    blah = vertcat(dataL1.ROI{i}{:});
%     blah = ceil(blah+trans);
    polyin = polyshape(blah(:,1),blah(:,2));
    [x,y] = centroid(polyin);
    ROIcentroidL1(i,:) = [x y];
    for ii = 1:size(blah,1)
        ROInew{i,1}{ii,1} = blah(ii,:);
    end
end
ROIcentroidL5 = [];ROInew = [];
for i = 1:length(dataL5.ROI)
    trans = rand(1)*350;
    blah = vertcat(dataL5.ROI{i}{:});
%     blah = ceil(blah+trans);
    polyin = polyshape(blah(:,1),blah(:,2));
    [x,y] = centroid(polyin);
    ROIcentroidL5(i,:) = [x y];
    for ii = 1:size(blah,1)
        ROInew{i,1}{ii,1} = blah(ii,:);
    end
end
%%
coorL5 = correlation_dice(Spikes);
Connected_ROI = Connectivity_dice(corr,0.3);
Connected_ROIL5 = Connectivity_dice(coorL5,0.1);
figure,Dendritic_Map(AverageImage,Connected_ROI,ROIcentroid,ROI,1,1),title('Pairwise Coopertivity')
figure,Dendritic_Map(AverageImage,Connected_ROIL5,ROIcentroidL5,dataL5.ROI,1,1),title('Pairwise Coopertivity')

%% 3d Plt
figure,
for i = 1:length(ROI)
    blah = vertcat(ROI{i}{:});
    a = 60; %PSF range
    b = 0; 
    adjustValue = ceil((b-a).*rand(1) + a);
    blah = [blah,adjustValue*ones(size(blah,1),1)];
    line(blah(:,1),blah(:,2),blah(:,3));
end
view(-30,20)

for i = 1:length(ROI)
    blah = vertcat(ROI{i}{:});
    a = -400; %PSF range
    b = -407; 
    adjustValue = ceil((b-a).*rand(1) + a);
    blah = [blah,adjustValue*ones(size(blah,1),1)];
    line(blah(:,1),blah(:,2),blah(:,3));
end
view(-30,20)













