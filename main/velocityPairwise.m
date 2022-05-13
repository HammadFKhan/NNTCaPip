function velocityPairwise(VR_data,Spikes)
%% Pairwise synchronization as velocity
Velocity = encoderVelocity(VR_data); % Bin width of 4 (ie. average every 4 readings)
% Extract window during velocity threshold
thresh = find(abs(Velocity(:,2))>1);
thresh = thresh*4; %multiply by bin width of the encoder time
% Now match threshold to the calcium time by window
for i = 1:length(thresh)
    CaTime = find(round(time)==thresh(i));
    syncWin{i} = Spikes(:,CaTime-30:CaTime+30);
    corr{i} = correlation_dice(syncWin{i});
    figure(i),htmp(corr{i});caxis([0.0 .5]);
end
%total pairwise Synchronization
velCorr = correlation_dice(horzcat(syncWin{:}));
figure,htmp(velCorr);caxis([0.0 .5]);
% Connectivity analysis during locomotion
%pairwise synchornization during rest
threshRest = find(Velocity(:,2)<1);
threshRest = threshRest*4; %multiply by bin width of the encoder time
% Now match threshold to the calcium time by window
for i = 1:length(threshRest)
    CaTime = find(round(time)==threshRest(i));
    syncWinRest{i} = Spikes(:,CaTime-30:CaTime+30);
end
%total pairwise Synchronization
velCorrRest = correlation_dice(horzcat(syncWinRest{:}));
figure,htmp(velCorrRest);caxis([0.0 .5]);
% Connectivity analyis at rest