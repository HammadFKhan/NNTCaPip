function [Ensemble] = ensembleStat(Ensemble)
% Ensemble stats for trial processing
for i = 1:size(Ensemble.NodeList,2)
    ensembleSize(i) = length(Ensemble.NodeList{i});
end
ensembleSize(ensembleSize==1) = [];
ensembleNum = length(ensembleSize);
ensembleRecruitment = (mean(ensembleSize)/max(cell2mat(Ensemble.NumNodes)));
minRecruitment = (min(ensembleSize)/max(cell2mat(Ensemble.NumNodes)));
maxRecruitment = (max(ensembleSize)/max(cell2mat(Ensemble.NumNodes)));


try % Stupid edge case when there is a string in cell (legacy noob build)
    activityCentroid = vertcat(Ensemble.ActivityCentroid{:});
    activityCentroidVariance = vertcat(Ensemble.ActivityCentroidVariance{:});
catch
    for i = 1:length(Ensemble.ActivityCentroid)
        cBuff{i} = class(Ensemble.ActivityCentroid{i}); % Find char class cells
        cvBuff{i} = class(Ensemble.ActivityCentroidVariance{i});
    end
    cIdx = find(~contains(cBuff,'char')); % Remove character cell arrays
    cvIdx = find(~contains(cvBuff,'char'));
    activityCentroid = vertcat(Ensemble.ActivityCentroid{cIdx}); % Do the analysis again
    activityCentroidVariance = vertcat(Ensemble.ActivityCentroidVariance{cvIdx});
end

% Output to structure
Ensemble.ensembleSize = ensembleSize;
Ensemble.ensembleNum = ensembleNum;
Ensemble.ensembleRecruitment = ensembleRecruitment;
Ensemble.minRecruitment = minRecruitment;
Ensemble.maxRecruitment = maxRecruitment;
Ensemble.activityCentroid = activityCentroid;
Ensemble.activityCentroidVariance = activityCentroidVariance;