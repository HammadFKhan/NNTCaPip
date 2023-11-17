% function ensembleVid(Ensemble,AverageImage,ROIcentroid,fn)
% if nargin<4 || strcmp('',fn)
%     fn = 'ensembleVideo';
% else
%     [~,name,~] = fileparts(fn)
%     fn = name;
% end
fn = 'ls_centroid'

set(0,'DefaultFigureWindowStyle','normal')
clear mov
writerobj = VideoWriter([fn '.avi'],'Uncompressed AVI'); % Initialize movie file
writerobj.FrameRate = 30;
open(writerobj);
[grad,~]=colorGradient([1 1 1] ,[1 0 0 ],55);
count = 1;
idx = find(~cellfun(@isempty,lateSpikeEnsemble.ActivityCentroid));
for j  = idx
    for i = 1:30
        fig = figure(count);title(['Late spike trial: ' num2str(j)]);
        axis([0 512 0 512])
        axis off
        set(gcf,'Position',[100 100 600 600])
        % Draw centroid area
        cx = lateSpikeEnsemble.ActivityCentroid{j}(1);
        cy = lateSpikeEnsemble.ActivityCentroid{j}(2);
        a = lateSpikeEnsemble.ActivityCentroidVariance{j}(1);
        b = lateSpikeEnsemble.ActivityCentroidVariance{j}(2);
        angle = 0;
        plot_ellipse(a,b,cx,cy,angle,grad(i,:)), hold on
        for ii = 1:length(ROI)
            blah = vertcat(ROI{ii}{:});
            line(blah(:,1),blah(:,2),'Color','black');
        end
        drawnow
        mov(count) = getframe(gcf);
        writeVideo(writerobj,mov(count));
        count = count+1;
        close(fig) 
    end
end
close(writerobj)
disp('Video saved to current directory')
