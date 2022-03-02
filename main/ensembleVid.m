function ensembleVid(Ensemble,AverageImage,ROIcentroid,fn)
if nargin<4 || strcmp('',fn)
    fn = 'ensembleVideo';
else
    [~,name,~] = fileparts(fn)
    fn = name;
end

set(0,'DefaultFigureWindowStyle','normal')
clear mov
writerobj = VideoWriter([fn '.avi'],'Uncompressed AVI'); % Initialize movie file
writerobj.FrameRate = 8;
open(writerobj);
for i = 1:Ensemble.ensembleIndentified 
    figure(i)
    axis off
    color = jet(Ensemble.ensembleIndentified);
    EnsembleMap(AverageImage,ROIcentroid,Ensemble.NodeList{i},8,color(i,:))
    set(gcf,'Position',[100 100 500 500])
    drawnow
    mov(i) = getframe(gcf);
    writeVideo(writerobj,mov(i));
end
close(writerobj)
disp('Video saved to current directory')
