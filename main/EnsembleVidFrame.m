% function ensembleVid(Ensemble,AverageImage,ROIcentroid,fn)
% if nargin<4 || strcmp('',fn)
%     fn = 'ensembleVideo';
% else
%     [~,name,~] = fileparts(fn)
%     fn = name;
% end
fn = 'centroid'

set(0,'DefaultFigureWindowStyle','normal')
clear mov
writerobj = VideoWriter([fn '.avi'],'Uncompressed AVI'); % Initialize movie file
writerobj.FrameRate = 8;
open(writerobj);
[grad,~]=colorGradient([1 1 1] ,[1 0 0 ],31);
ensemblePlot = vertcat(lateSpikeEnsemble1.rankEnsembles{1:2},lateSpikeEnsemble2.rankEnsembles{1:2}...
    ,lateSpikeEnsemble3.rankEnsembles{1:2},lateSpikeEnsemble4.rankEnsembles{1:2});
ROIcentroid = vertcat(lateSpikeEnsemble1.ROIcentroid,lateSpikeEnsemble2.ROIcentroid,...
    lateSpikeEnsemble3.ROIcentroid,lateSpikeEnsemble4.ROIcentroid);
count = 1;
for i = 1:30
    figure(count)
    axis off
    EnsembleMap(AverageImage,ROIcentroid,ensemblePlot,6,grad(i,:))
    set(gcf,'Position',[100 100 600 600])
    drawnow
    mov(count) = getframe(gcf);
    writeVideo(writerobj,mov(count));
    count = count+1;
end
% EnsembleMap(AverageImage,ROIcentroid,ensemblePlot,6,[1 1 1])
[grad,~]=colorGradient([1 1 1] ,[0 0 1 ],31);
ensemblePlot = [];
ensemblePlot = vertcat(nolateSpikeEnsemble1.rankEnsembles{1:2},nolateSpikeEnsemble2.rankEnsembles{1:2}...
    ,nolateSpikeEnsemble3.rankEnsembles{1:2},nolateSpikeEnsemble4.rankEnsembles{1:2});
for i = 1:30
    figure(count)
    axis off
    EnsembleMap(AverageImage,ROIcentroid,ensemblePlot,6,grad(i,:))
    set(gcf,'Position',[100 100 600 600])
    drawnow
    mov(count) = getframe(gcf);
    writeVideo(writerobj,mov(count));
    count = count+1;
end
close(writerobj)
disp('Video saved to current directory')
