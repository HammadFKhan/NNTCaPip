%% Entropy Stats
%Entropy parser
clear
pathname = uigetdir(pwd,'Input Directory');
pathname = fullfile(pathname);
directory = dir(fullfile(pathname,'*.mat'));
L = length(directory);
for i = 1:L
    S = load(fullfile(directory(i).folder,directory(i).name));
    StrName=fieldnames(S);
    StrName = StrName{1};
    [P_c{i},I_c{i},H_c(i,1),H_e(i,1),sizeE{i},sizeEdge{i}] = EntropyParser(S.(StrName));
end 
%%
figure,hold on
x = 1.1*ones(length(H_c),1);
boxplot(H_c,'plotstyle','compact'), title('Entropy by Ensemble cell')
scatter(x,H_c,'filled','k')
xlim([0 2])
ylim([3 9])
figure,hold on
boxplot(H_e,'plotstyle','compact'), title('Entropy by Ensemble Frequency')
scatter(x,H_e,'filled','k')
xlim([0 2])
ylim([3 9])

%%
trialSize = cellfun(@size,sizeE,'UniformOutput',false)';
trialSize = cell2mat(trialSize);
trialSize(:,1) = [];
trialMax = max(trialSize);
figure
for i = 1:L
    if length(sizeE{i})<trialMax
        sizeE{i} = horzcat(sizeE{i},zeros(1,trialMax-length(sizeE{i})));
    end
    loglog(sizeE{i}),title('Ranked Ensembles'); hold on
end
errorE = vertcat(sizeE{:});
x1 = 1:size(errorE,2);
figure,lineError(x1,errorE,'ste');
xlim([0 100])
ylim([0 1])