% LFP + calcium analysis

%% Beta calcium coupling

startFrame = ceil(LFP.beta.betaBurst.detectedBeta(:,1)*CaFR);
endFrame = ceil(LFP.beta.betaBurst.detectedBeta(:,3)*CaFR);
peakFrame = ceil(LFP.beta.betaBurst.detectedBeta(:,2)*CaFR);
betaEventFrame = [startFrame endFrame peakFrame]; %calculate the calcium frame in which beta event window occurs

 % Calculate coactivity during beta event frame
 for i = 1:size(betaEventFrame,1)
     betaCaCoupling(i) = median(coactive_cells(:,betaEventFrame(i,1):betaEventFrame(i,2)));
 end
 
 % Calculate coactivity to beta event power

figure,plot(betaCaCoupling,LFP.beta.betaBurst.detectedBeta(:,4),'k.')
 
figure,histogram(betaCaCoupling,0:0.01:.1);box off

for i = 1:341
    xline(betaEventFrame(i,3),'r')
end
