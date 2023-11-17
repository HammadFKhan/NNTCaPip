function output = LFPmodulationIndex(DeltaFoverF,Ensemble,LFP)
%% Modulation probability
% we reasoned that limitations in calcium sensitivity, baseline levels and
% ripple amplitude can influence the response activity seen in different
% subtypes. To circumvent those limitations, we developed a measure of
% modulation to evaluate whether a cell had an increased activity during a
% SWR event statistically higher than its baseline fluctuation during
% immobility. To do so, we computed for each cell a shuffle distribution
% that consisted of 1000 values. Each value was the average of N randomly
% chosen ΔF/F points during immobility (N corresponds to the number of SWRs
% detected in the session). Then, for each SWR, if the difference in
% activity in a 400ms window centered around the onset was greater than the
% 95th percentile of the shuffle distribution, the cell was modulated
% during this event. The modulation probability reports the percentage of
% SWRs leading to an increase above baseline.
% Taken from https://doi.org/10.1016/j.neuron.2020.09.013

for i = 1:size(DeltaFoverF,1)
    for ii = 1:size(LFP.beta.betaBurst.detectedBeta,1)
        win = [floor(LFP.beta.betaBurst.detectedBeta(ii,2)*10.57)-20,floor(LFP.beta.betaBurst.detectedBeta(ii,2)*10.57)+20];
        try
        temp1 = DeltaFoverF(i,win(1):win(1)+10);
        temp2 = DeltaFoverF(i,win(2)-10:win(2));
        if sum(temp2-temp1)>0
            modulationIndex(ii) = 1;
        else
            modulationIndex(ii) = 0;
        end
        catch
            continue;
        end
    end
    modulationProbability(i) = mean(modulationIndex);
end
%% Now calculate it for each ensemble neuron in the ensemble

for i = 1:length(Ensemble.rankEnsembles)
    idx = Ensemble.rankEnsembles{i};
   ensembleModulationProbability{i} = modulationProbability(idx)';
end
combinedEnsemble = vertcat(ensembleModulationProbability{:});
ensembleModulation = cellfun(@(x) nanmean(x), ensembleModulationProbability);

pooledEnsembleModulationProbability = [];
for i = 1:length(ensembleModulationProbability)
    temp = [ensembleModulationProbability{i};zeros(length(ensembleModulationProbability{1})-length(ensembleModulationProbability{i}),1)];
    pooledEnsembleModulationProbability(:,i) = temp;
end
%% Correlation of ensemble size
shuf = randperm(24,24);
ensembleSize = cellfun(@(x) length(x), Ensemble.rankEnsembles);
[p,S] = polyfit(ensembleSize,ensembleModulation,1); 
y = ensembleModulation; 
f = polyval(p,ensembleSize); 
SStot = sum((y-mean(y)).^2);                    % Total Sum-Of-Squares
SSres = sum((y-f).^2);                       % Residual Sum-Of-Squares
Rsq = 1-SSres/SStot;
%% Outputs
output.modulationIndex = modulationIndex;
output.modulationProbability = modulationProbability;
output.ensembleModulationProbability = ensembleModulationProbability;
output.pooledEnsembleModulationProbability = pooledEnsembleModulationProbability;
output.p = p;
output.S = S;
output.Rsq = Rsq;