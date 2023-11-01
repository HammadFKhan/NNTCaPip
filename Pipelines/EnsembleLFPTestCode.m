%% Test code for L23 -> L5 ensemble mapping during beta/gamma modulation
ROIcentroid = [];
for i = 1:length(ROI)
    blah = vertcat(ROI{i}{:});
    ROIcentroid(i,:) = floor(mean(blah,1));
end

%% Spike detection from dF/F
std_threshold = 3;      % from Carrilo-Reid and Jordan Hamm's papers
static_threshold = .01;
L5.Spikes = rasterizeDFoF(L5.DeltaFoverF,std_threshold,static_threshold);
L23.Spikes = rasterizeDFoF(L23.DeltaFoverF,std_threshold,static_threshold);
%% Ensemble Analysis
factorCorrection = 5*floor(size(L23.Spikes,2)/5); % Correct for frame size aquisition
L23.Ensemble = ensembleAnalysis(L23.Spikes(:,1:factorCorrection),L23.ROIcentroid);
L23.Ensemble = ensembleMetric(L23.Ensemble,L23.AverageImage,L23.ROIcentroid);
L23.Ensemble = ensembleStat(L23.Ensemble);

L5.Ensemble = ensembleAnalysis(L5.Spikes(:,1:factorCorrection),L5.ROIcentroid);
L5.Ensemble = ensembleMetric(L5.Ensemble,L5.AverageImage,L5.ROIcentroid);
L5.Ensemble = ensembleStat(L5.Ensemble);
%% Stability Analysis
L23ensembleStability = vertcat(L23.Ensemble.ensembleStability{:});
L23ensembleStability(L23ensembleStability==0) = NaN;
L5ensembleStability = vertcat(L5.Ensemble.ensembleStability{:});
L5ensembleStability(L5ensembleStability==0) = NaN;
%% Centroid Analysis
L23Centroid = sqrt(L23.Ensemble.activityCentroidVariance(:,1).^2+L23.Ensemble.activityCentroidVariance(:,2).^2);
L5Centroid = sqrt(L5.Ensemble.activityCentroidVariance(:,1).^2+L5.Ensemble.activityCentroidVariance(:,2).^2);
%% Ensemble Activation Tracking across Layers
% L23->L5
[L23vectorized,~] = cosine_similarity(L23.Spikes(L23.Ensemble.rankEnsembles{1},1:factorCorrection),5);
[L5vectorized,~] = cosine_similarity(L5.Spikes(L5.Ensemble.rankEnsembles{1},1:factorCorrection),5);
%% Build data for ensemble causality
ensembleNumInfer = 20; % Number of ensembles to infer from for each state
l23Caus = [];l5Caus = [];
for n = 1:ensembleNumInfer
l23Caus = vertcat(l23Caus,mean(L23.DeltaFoverF(L23.Ensemble.rankEnsembles{n},:)));
l5Caus = vertcat(l5Caus,mean(L5.DeltaFoverF(L5.Ensemble.rankEnsembles{n},:)));
end
dat = [l23Caus;l5Caus];
ensembleCausality = mvgcstate(dat);
%% Causality stats
% Now find the causality flow direction for each pair significant values
% this is a bit tricky since we'll have to reference the original dat
% build. In this case we did 20/20 (ensembleNumInfer) L23 L5 ensemble mapping so we use these
% as indexes. We'll include the probability for inferred causal
% relationship and the associated GC coefficient (for additional Ns)
[row,col] = find(ensembleCausality.sig==1);
sigIdx = [row col]; % Only grab data points from significant relationships
L23L23F = [];L23L5F = [];L5L5F = [];L5L23F = [];
L23L23Coor = [];L23L5Coor = [];L5L5Coor = [];L5L23Coor = [];
for n = 1:size(sigIdx,1)
    if sigIdx(n,1)<=ensembleNumInfer && sigIdx(n,2)<=ensembleNumInfer
        L23L23Coor = vertcat(L23L23Coor,sigIdx(n,:));
        L23L23F = [L23L23F;ensembleCausality.Fnew(sigIdx(n,1),sigIdx(n,2))];
    elseif sigIdx(n,1)<=ensembleNumInfer && sigIdx(n,2)>ensembleNumInfer
        L23L5Coor = vertcat(L23L5Coor,sigIdx(n,:));
        L23L5F = [L23L5F;ensembleCausality.Fnew(sigIdx(n,1),sigIdx(n,2))];
    elseif sigIdx(n,1)>ensembleNumInfer && sigIdx(n,2)>ensembleNumInfer
        L5L5Coor = vertcat(L5L5Coor,sigIdx(n,:));
        L5L5F = [L5L5F;ensembleCausality.Fnew(sigIdx(n,1),sigIdx(n,2))];
    else
        L5L23Coor = vertcat(L5L23Coor,sigIdx(n,:));
        L5L23F = [L5L23F;ensembleCausality.Fnew(sigIdx(n,1),sigIdx(n,2))];
    end
end
ensembleCausality.stats.GC.L23L23 = L23L23F;
ensembleCausality.stats.GC.L23L5 = L23L5F;
ensembleCausality.stats.GC.L5L5 = L5L5F;
ensembleCausality.stats.GC.L5L23 = L5L23F;
ensembleCausality.stats.L23L23Coor = L23L23Coor;
ensembleCausality.stats.L23L5Coor = L23L5Coor;
ensembleCausality.stats.L5L5Coor = L5L5Coor;
ensembleCausality.stats.L5L23Coor = L5L23Coor;
%% Ensemble structure during causal (inferred) activation states
% We'll be doing this by calculating the structure of ensembles 
% based on their inferred activation direction. For instance L23-> L23
% activation of all ensemble pairs and their coinciding ensemble stability
% metric. This will give us an overview of volumetric cortical activation
% flow.
% I made a little function on the bottom for redundancy

% Make ensemble stability array for conditions
padding = size(L23ensembleStability,2)-size(L5ensembleStability,2);
if padding<0
    ensembleCausality.causalEnsembleStability =...
        vertcat(padarray(L23ensembleStability(1:ensembleNumInfer,:),[1 abs(padding)],NaN,'post'),L5ensembleStability(1:ensembleNumInfer,:));
else
    ensembleCausality.causalEnsembleStability =...
        vertcat(L23ensembleStability(1:ensembleNumInfer,:),padarray(L5ensembleStability(1:ensembleNumInfer,:),[1 abs(padding)],NaN,'post'));
end

% L23->L23
ensembleCausality.stats.ensembleStability.L23L23 =...
    causalStability(ensembleCausality.stats.L23L23Coor,ensembleCausality.causalEnsembleStability);
% L23->L5 
ensembleCausality.stats.ensembleStability.L23L5 = ...
    causalStability(ensembleCausality.stats.L23L5Coor,ensembleCausality.causalEnsembleStability);
% L5->L5
ensembleCausality.stats.ensembleStability.L5L5 = ...
    causalStability(ensembleCausality.stats.L5L5Coor,ensembleCausality.causalEnsembleStability);
% L5->L23
ensembleCausality.stats.ensembleStability.L5L23 = ...
    causalStability(ensembleCausality.stats.L5L23Coor,ensembleCausality.causalEnsembleStability);
%% Plot out ensemble stability during inferred causal conditions
figure,customBoxplot(ensembleCausality.stats.ensembleStability.L23L23),box off, set(gca,'TickDir','out'),set(gca,'FontSize',16),title('L23->L23')
figure,customBoxplot(ensembleCausality.stats.ensembleStability.L23L5),box off, set(gca,'TickDir','out'),set(gca,'FontSize',16),title('L23->L5')
figure,customBoxplot(ensembleCausality.stats.ensembleStability.L5L5),box off, set(gca,'TickDir','out'),set(gca,'FontSize',16),title('L5->L5')
figure,customBoxplot(ensembleCausality.stats.ensembleStability.L5L23),box off, set(gca,'TickDir','out'),set(gca,'FontSize',16),title('L5->L23')
%% Now start looking at LFP relationships
% For now we'll look at LFP power envelopes across narrow bands during
% periods of ensemble activations. Which we can extract from population

% Calculate narrowband LFP envelopes
% coactivity of ensembles across the recording period
LFP.wideBandEnvelope = calcEnvelope(LFP.dsWideband);
LFP.betaEnvelope = calcEnvelope(LFP.dsbetaLFP);
LFP.gammaEnvelope = calcEnvelope(LFP.dsgammaLFP);

figure,subplot(4,1,1),plot(time,smoothdata(sum(L23.Spikes(L23.Ensemble.rankEnsembles{1},:)),'gaussian',30));
subplot(4,1,2),plot(LFP.dsTime,LFP.wideBandEnvelope),hold on,plot(LFP.dsTime,LFP.dsWideband)
subplot(4,1,3),plot(LFP.dsTime,LFP.betaEnvelope),hold on,plot(LFP.dsTime,LFP.dsbetaLFP)
subplot(4,1,4),plot(LFP.dsTime,LFP.gammaEnvelope),hold on,plot(LFP.dsTime,LFP.dsgammaLFP)
%% find coactive periods of ensemble activation
% Built window response TODO: window response needs to be sweeped for
% smoothing effects and how it effects ensemble activation
L23 = ensembleNarrowbandCoupling(L23,LFP);
L5 = ensembleNarrowbandCoupling(L5,LFP);
figure,customBoxplot([L23.LFP.betaEnsembleCoefficients',L23.LFP.gammaEnsembleCoefficients']),title('L23')
figure,customBoxplot([L5.LFP.betaEnsembleCoefficients',L5.LFP.gammaEnsembleCoefficients']),title('L5')
%% local functions
function output = causalStability(ensembleCausality,ensembleStability)
for n = 1:size(ensembleCausality,1)
    output(n,:) = nanmean(ensembleStability(ensembleCausality(n,:),:),2)';
end
end

function output = calcEnvelope(lfp)
[upper,lower] = envelope(lfp,100,'peak');
output = abs(upper)+abs(lower);
end