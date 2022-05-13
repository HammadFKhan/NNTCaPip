clear
pathname = uigetdir(pwd,'Input Directory');
pathname = fullfile(pathname);
directory = dir(fullfile(pathname,'*.mat'));
L = length(directory);
for i = 1:L
    disp(['Parsing: ' num2str(directory(i).name)])
    S = load(fullfile(directory(i).folder,directory(i).name));
    StrName=fieldnames(S);
    StrName = StrName{1};
%     [P_c{i},I_c{i},H_c(i,1),H_e(i,1),sizeE{i},sizeEdge{i}] = EntropyParser(S.(StrName));
    [ensembleSize{i},ensembleNum{i},ensembleRecruitment{i},minRecruitment{i},maxRecruitment{i}] = ensembleStat(S.(StrName));
end 
%%
figure,bar(cell2mat(ensembleRecruitment))
figure,bar(cell2mat(maxRecruitment))
figure,bar(cell2mat(minRecruitment))
figure,bar(cell2mat(ensembleNum))
for i = 1:33
    ensembleSizeMean(i) = mean(ensembleSize{i});
end
%%
data6 = test(25:end);
figure,boxplot(data,'PlotStyle','compact'),ylim([0 0.5]),box off,hold on
plot(1.1*ones(1,length(data)),data,'.');