function [ensembleSize,ensembleNum,ensembleRecruitment,minRecruitment,maxRecruitment] = ensembleStat(Ensemble,ROI)
% Ensemble stats for trial processing
for i = 1:size(Ensemble.NodeList,2)
    ensembleSize(i) = length(Ensemble.NodeList{i});
end
ensembleSize(ensembleSize==1) = [];
ensembleNum = length(ensembleSize);
ensembleRecruitment = (mean(ensembleSize)/size(ROI,1));
minRecruitment = (min(ensembleSize)/size(ROI,1));
maxRecruitment = (max(ensembleSize)/size(ROI,1));

