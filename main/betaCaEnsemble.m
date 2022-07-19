function stateLFP = betaCaEnsemble(Spikes,SpikesFrame,Ensemble,LFP,CaFR)

% Input variables: 
% Spikes: NxT array to calculate the coactivity index of the cells
% Spikes Frame: the correct frame in which the ensemble events occured
% LFP: structure array for single channel LFP data
% CaFR: frame rate of imaging aquisition

coactive_cells = mean(Spikes,1);

startFrame = ceil(LFP.beta.betaBurst.detectedBeta(:,1)*CaFR);
endFrame = ceil(LFP.beta.betaBurst.detectedBeta(:,3)*CaFR);
peakFrame = ceil(LFP.beta.betaBurst.detectedBeta(:,2)*CaFR);
betaEventFrame = [startFrame endFrame peakFrame]; %calculate the calcium frame in which beta event window occurs

% Beta events within ensembles
frameT = SpikesFrame(Ensemble.ensembleFrame); %set the correct index from ensemble->stateSpikes->Spikes
count = 1;
for i = 1:size(Ensemble.ensemble,2) %state dependant ensemble
    temp = find(frameT(i)==betaEventFrame(:,3)); % checks to see if a beta event lies on the frame\
    if temp
        disp(['Beta event matched to idx: ' num2str(i)])
        hold on,xline(i,'r');
        ensembleBetaMatch(count,1) = frameT(i); %frame
        ensembleBetaMatch(count,2) = temp; %beta index
        count = count+1;
    end
end

% Calculate coactivity during beta event frame
betaCaCoupling = coactive_cells(:,ensembleBetaMatch(:,2));

stateLFP.beta = LFP.beta; % create a structure array looking at only behavior based beta
stateLFP.beta.ensembleBetaMatch = ensembleBetaMatch;
stateLFP.beta.betaBurst.detectedBeta = LFP.beta.betaBurst.detectedBeta(ensembleBetaMatch(:,2),:);
stateLFP.beta.betaCaCoupling = betaCaCoupling;