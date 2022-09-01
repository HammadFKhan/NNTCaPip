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
[grad,~]=colorGradient([0 0 1] ,[0 0 0],Ensemble.ensembleIndentified+1)
for i = 1:30
    figure(i)
    axis off
    EnsembleMap(AverageImage,ROIcentroid,Ensemble.rankEnsembles{i},6,grad(i,:))
    set(gcf,'Position',[100 100 600 600])
    drawnow
    mov(i) = getframe(gcf);
    writeVideo(writerobj,mov(i));
end
close(writerobj)
disp('Video saved to current directory')
