function Calcium = ensembleNarrowbandCoupling(Calcium,LFP)
% Calcium data structure needs to have Spikes and Ensemble structure mapped
% from previous functions
% LFP needs to have narrowband envelope response
% Built window response TODO: window response needs to be sweeped for
% smoothing effects and how it effects ensemble activation
ensembleNum = Calcium.Ensemble.ensembleNum;
gammaEnsemble = cell(1,ensembleNum);betaEnsemble = cell(1,ensembleNum);
fprintf('Calculating ensemble-LFP coupling...')
for nn = 1:ensembleNum
ensembleAct = sum(Calcium.Spikes(Calcium.Ensemble.rankEnsembles{nn},:));
idxE = find(ensembleAct>1); % divide by calcium sample rate for locking with LFP
ensembleAct = smoothdata(ensembleAct(idxE),'gaussian',10);
idxE = idxE/10.5759; %Ca Fs

winlen = .4; % arbitrary chosen
for n = 1:length(idxE)
    win = ceil((idxE(n)-winlen)*LFP.downsampledFs:(idxE(n)+winlen)*LFP.downsampledFs);
    % Discrete ensemble activation value we then integrate under LFP
    % envelope
    gammaEnsemble{nn} = vertcat(gammaEnsemble{nn},mean(LFP.gammaEnvelope(win)));
    betaEnsemble{nn} = vertcat(betaEnsemble{nn}, mean(LFP.betaEnvelope(win)));
end
% Build coefficient linear model of data to find ensembles that are
% correlated to beta activity
mdl = fitlm(ensembleAct,betaEnsemble{nn});
betaEnsembleCoefficients(nn) = mdl.Rsquared.Ordinary;
mdl = fitlm(ensembleAct,gammaEnsemble{nn});
gammaEnsembleCoefficients(nn) = mdl.Rsquared.Ordinary;
end
fprintf('done\n')
Calcium.LFP.betaEnsembleCoefficientCutoff = 1.5*nanstd(betaEnsembleCoefficients)+nanmean(betaEnsembleCoefficients);
Calcium.LFP.gammaEnsembleCoefficientCutoff = 1.5*nanstd(gammaEnsembleCoefficients)+nanmean(gammaEnsembleCoefficients);
Calcium.LFP.betaEnsembleCoefficientsRaw = betaEnsembleCoefficients;
Calcium.LFP.gammaEnsembleCoefficientsRaw = gammaEnsembleCoefficients;
betaEnsembleCoefficients(betaEnsembleCoefficients<Calcium.LFP.betaEnsembleCoefficientCutoff) = NaN;
gammaEnsembleCoefficients(gammaEnsembleCoefficients<Calcium.LFP.gammaEnsembleCoefficientCutoff) = NaN;

Calcium.LFP.betaEnsembleCoefficients = betaEnsembleCoefficients;
Calcium.LFP.gammaEnsembleCoefficients = gammaEnsembleCoefficients;