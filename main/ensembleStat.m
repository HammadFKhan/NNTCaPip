function [ensembleSize,ensembleNum,ensembleRecruitment,minRecruitment,maxRecruitment] = ensembleStat(Ensemble)
% Ensemble stats for trial processing
for i = 1:size(Ensemble.NodeList,2)
    ensembleSize(i) = length(Ensemble.NodeList{i});
end
ensembleSize(ensembleSize==1) = [];
ensembleNum = length(ensembleSize);
ensembleRecruitment = (mean(ensembleSize)/max(cell2mat(Ensemble.NumNodes)));
minRecruitment = (min(ensembleSize)/max(cell2mat(Ensemble.NumNodes)));
maxRecruitment = (max(ensembleSize)/max(cell2mat(Ensemble.NumNodes)));

