function DeltaFoverFplotter(DeltaFoverF,stdev_threshold,static_threshold)


for i = 1:size(DeltaFoverF,3)
    num_images = size(DeltaFoverF,2)+1;
    subplot(1,size(DeltaFoverF,3),i)
    DeltaFplot = DeltaFoverF(:,:,i);
    plot(linspace(1,num_images-1,num_images-1),DeltaFplot)
    hold on
    plot(linspace(1,num_images-1,num_images-1),(zeros(1,num_images-1)+stdev_threshold*std(std(DeltaFoverF(:,:,i)))));
    plot(linspace(1,num_images-1,num_images-1),(zeros(1,num_images-1)+static_threshold));
    ylim ([-1.1*max(max(max(DeltaFoverF))),1.1*max(max(max(DeltaFoverF)))])
end