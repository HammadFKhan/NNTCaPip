function Calcium = leverCaModulation(DeltaFoverF,Spikes,Behaviour,caTime)
%Genenerate behavior triggered calcium response
if ~exist('caTime','var')
    error('Include calcium time using volume frame rate provided through the metaData!')
end
% Make cutoff for array sake (may drop one frame in some trials which is
% acceptable for us)
caTimeCutoff = floor((Behaviour.parameters.windowBeforeCue*diff(caTime(1:2))^-1+Behaviour.parameters.windowAfterCue*diff(caTime(1:2))^-1));

Calcium = struct();
caM = mean(DeltaFoverF,2);
for i = 1:length(Behaviour.cueHitTrace)
    %Calculate indices of interest from calcium data; can do this using Intan time and caFR
    caIndex = find(caTime>=floor(Behaviour.cueHitTrace(i).LFPtime(1)) & caTime<=floor(Behaviour.cueHitTrace(i).LFPtime(end)));
    if length(caIndex)>caTimeCutoff
        caIndex = caIndex(1:caTimeCutoff);
    end
    Calcium.hit.caTrace(i).DeltaFoverF = DeltaFoverF(:,caIndex);
    Calcium.hit.modulationIdx(:,i) = single(1.2*mean(Calcium.hit.caTrace(i).DeltaFoverF,2)>caM);
    Calcium.hit.modulationIdxM = mean(Calcium.hit.modulationIdx,2);
    Calcium.hit.DeltaFoverF(:,:,i) = Calcium.hit.caTrace(i).DeltaFoverF;
    Calcium.hit.Spikes(:,:,i) = Spikes(:,caIndex);
end
for i = 1:length(Behaviour.cueMissTrace)
    caIndex = find(caTime>=floor(Behaviour.cueMissTrace(i).LFPtime(1)) & caTime<=floor(Behaviour.cueMissTrace(i).LFPtime(end)));
    if length(caIndex)>caTimeCutoff
        caIndex = caIndex(1:caTimeCutoff);
    end
    Calcium.miss.caTrace(i).DeltaFoverF = DeltaFoverF(:,caIndex);
    Calcium.miss.modulationIdx(:,i) = single(1.2*mean(Calcium.miss.caTrace(i).DeltaFoverF,2)>caM);
    Calcium.miss.modulationIdxM = mean(Calcium.miss.modulationIdx,2);
    Calcium.miss.DeltaFoverF(:,:,i) = Calcium.miss.caTrace(i).DeltaFoverF;
    Calcium.miss.Spikes(:,:,i) = Spikes(:,caIndex);
end
if isfield(Behaviour,'MIHitTrace')
    for i = 1:length(Behaviour.MIHitTrace)
        caIndex = find(caTime>=floor(Behaviour.MIHitTrace(i).LFPtime(1)) & caTime<=floor(Behaviour.MIHitTrace(i).LFPtime(end)));
        if length(caIndex)>caTimeCutoff
            caIndex = caIndex(1:caTimeCutoff);
        end
        Calcium.MIhit.caTrace(i).DeltaFoverF = DeltaFoverF(:,caIndex);
        Calcium.MIhit.modulationIdx(:,i) = single(1.2*mean(Calcium.MIhit.caTrace(i).DeltaFoverF,2)>caM);
        Calcium.MIhit.modulationIdxM = mean(Calcium.MIhit.modulationIdx,2);
        Calcium.MIhit.DeltaFoverF(:,:,i) = Calcium.MIhit.caTrace(i).DeltaFoverF;
        Calcium.MIhit.Spikes(:,:,i) = Spikes(:,caIndex);
    end
end
if isfield(Behaviour,'MIFATrace')
    for i = 1:length(Behaviour.MIFATrace)
        caIndex = find(caTime>=floor(Behaviour.MIFATrace(i).LFPtime(1)) & caTime<=floor(Behaviour.MIFATrace(i).LFPtime(end)));
        if length(caIndex)>caTimeCutoff
            caIndex = caIndex(1:caTimeCutoff);
        end
        Calcium.MIFA.caTrace(i).DeltaFoverF = DeltaFoverF(:,caIndex);
        Calcium.MIFA.modulationIdx(:,i) = single(1.2*mean(Calcium.MIFA.caTrace(i).DeltaFoverF,2)>caM);
        Calcium.MIFA.modulationIdxM = mean(Calcium.MIFA.modulationIdx,2);
        Calcium.MIFA.DeltaFoverF(:,:,i) = Calcium.MIFA.caTrace(i).DeltaFoverF;
        Calcium.MIFA.Spikes(:,:,i) = Spikes(:,caIndex);
    end
end

% Define motor neurons (modulation during either hit/miss trials)
Calcium.motorNeuron = unique(find(Calcium.hit.modulationIdxM>0.5 | Calcium.miss.modulationIdxM>0.5 |Calcium.MIFA.modulationIdxM>0.5));
% Define hit/reward neurons (modulated stronger during hit trials)
Calcium.hitNeuron = unique(find(Calcium.hit.modulationIdxM>0.5 & Calcium.hit.modulationIdxM>Calcium.miss.modulationIdxM));
% Define miss neurons (modulated stronger during miss trials)
Calcium.missNeuron = unique(find(Calcium.miss.modulationIdxM>0.5 & Calcium.miss.modulationIdxM>Calcium.hit.modulationIdxM));
% Define FA neurons (modulated stronger during FA trials (absence of cue) )
Calcium.FANeuron = unique(find(Calcium.MIFA.modulationIdxM>0.5 & Calcium.MIFA.modulationIdxM>Calcium.hit.modulationIdxM));

% Build lever task calcium traces
for i = 1:length(Behaviour.cueHitTrace)
    % Do for motor neuron
    [~,caMaxIdx] = max(Calcium.hit.caTrace(i).DeltaFoverF(:,:),[],2);
    [~,idx] = sort(caMaxIdx);
    motorNeuronTrace = Calcium.hit.caTrace(i).DeltaFoverF(:,:);
    Calcium.hit.motorNeuronTrace(:,:,i) = motorNeuronTrace(idx,:);
    
    % Do for hit neuron
    [~,caMaxIdx] = max(Calcium.hit.caTrace(i).DeltaFoverF(:,:),[],2);
    [~,idx] = sort(caMaxIdx);
    hitNeuronTrace = Calcium.hit.caTrace(i).DeltaFoverF(:,:);
    Calcium.hit.hitNeuronTrace(:,:,i) = hitNeuronTrace(idx,:);
    
    % Do for miss neuron
    [~,caMaxIdx] = max(Calcium.hit.caTrace(i).DeltaFoverF(Calcium.missNeuron,:),[],2);
    [~,idx] = sort(caMaxIdx);
    missNeuronTrace = Calcium.hit.caTrace(i).DeltaFoverF(Calcium.missNeuron,:);
    Calcium.hit.missNeuronTrace(:,:,i) = missNeuronTrace(idx,:);
end
for i = 1:length(Behaviour.cueMissTrace)
    % Do for motor neuron
    [~,caMaxIdx] = max(Calcium.miss.caTrace(i).DeltaFoverF(Calcium.motorNeuron,:),[],2);
    [~,idx] = sort(caMaxIdx);
    motorNeuronTrace = Calcium.miss.caTrace(i).DeltaFoverF(Calcium.motorNeuron,:);
    Calcium.miss.motorNeuronTrace(:,:,i) = motorNeuronTrace(idx,:);
    
    % Do for hit neuron
    [~,caMaxIdx] = max(Calcium.miss.caTrace(i).DeltaFoverF(Calcium.hitNeuron,:),[],2);
    [~,idx] = sort(caMaxIdx);
    hitNeuronTrace = Calcium.miss.caTrace(i).DeltaFoverF(Calcium.hitNeuron,:);
    Calcium.miss.hitNeuronTrace(:,:,i) = hitNeuronTrace(idx,:);
    
    % Do for miss neuron
    [~,caMaxIdx] = max(Calcium.miss.caTrace(i).DeltaFoverF(Calcium.missNeuron,:),[],2);
    [~,idx] = sort(caMaxIdx);
    missNeuronTrace = Calcium.miss.caTrace(i).DeltaFoverF(Calcium.missNeuron,:);
    Calcium.miss.missNeuronTrace(:,:,i) = missNeuronTrace(idx,:);
end
for i = 1:length(Behaviour.MIFATrace)
    % Do for motor neuron
    [~,caMaxIdx] = max(Calcium.MIFA.caTrace(i).DeltaFoverF(Calcium.motorNeuron,:),[],2);
    [~,idx] = sort(caMaxIdx);
    motorNeuronTrace = Calcium.MIFA.caTrace(i).DeltaFoverF(Calcium.motorNeuron,:);
    Calcium.MIFA.motorNeuronTrace(:,:,i) = motorNeuronTrace(idx,:);
    
    % Do for hit neuron
    [~,caMaxIdx] = max(Calcium.MIFA.caTrace(i).DeltaFoverF(Calcium.hitNeuron,:),[],2);
    [~,idx] = sort(caMaxIdx);
    hitNeuronTrace = Calcium.MIFA.caTrace(i).DeltaFoverF(Calcium.hitNeuron,:);
    Calcium.MIFA.hitNeuronTrace(:,:,i) = hitNeuronTrace(idx,:);
    
    % Do for miss neuron
    [~,caMaxIdx] = max(Calcium.MIFA.caTrace(i).DeltaFoverF(Calcium.missNeuron,:),[],2);
    [~,idx] = sort(caMaxIdx);
    missNeuronTrace = Calcium.MIFA.caTrace(i).DeltaFoverF(Calcium.missNeuron,:);
    Calcium.MIFA.missNeuronTrace(:,:,i) = missNeuronTrace(idx,:);
end

%% plot some stuff
plotOn = 1;
if plotOn
    %figure,boxplot([Calcium.hit.modulationIdxM,Calcium.miss.modulationIdxM])
%     figure,plot(ones(1,size(caM,1)),Calcium.hit.modulationIdxM,'k.'),hold on,title('All modulated neurons')
%     plot(2*ones(1,size(caM,1)),Calcium.miss.modulationIdxM,'k.'),hold on
%     for i = 1:size(caM,1)
%         line([1 2],[Calcium.hit.modulationIdxM(i) Calcium.miss.modulationIdxM(i)]),xlim([0.5 2.5])
%     end
%     
%     figure,plot(ones(1,size(Calcium.motorNeuron,1)),Calcium.hit.modulationIdxM(Calcium.motorNeuron),'k.')...
%         ,hold on,title('Motor Neuron')
%     plot(2*ones(1,size(Calcium.motorNeuron,1)),Calcium.miss.modulationIdxM(Calcium.motorNeuron),'k.'),hold on
%     for i = Calcium.motorNeuron
%         line([1 2],[Calcium.hit.modulationIdxM(i) Calcium.miss.modulationIdxM(i)]),xlim([0.5 2.5])
%     end
%     
%     figure,plot(ones(1,size(Calcium.hitNeuron,1)),Calcium.hit.modulationIdxM(Calcium.hitNeuron),'k.')...
%         ,hold on,title('Hit Responsive Neuron')
%     plot(2*ones(1,size(Calcium.hitNeuron,1)),Calcium.miss.modulationIdxM(Calcium.hitNeuron),'k.'),hold on
%     for i = Calcium.hitNeuron
%         line([1 2],[Calcium.hit.modulationIdxM(i) Calcium.miss.modulationIdxM(i)]),xlim([0.5 2.5])
%     end
%     
%     figure,plot(ones(1,size(Calcium.missNeuron,1)),Calcium.hit.modulationIdxM(Calcium.missNeuron),'k.')...
%         ,hold on,title('Miss Responsive Neuron')
%     plot(2*ones(1,size(Calcium.missNeuron,1)),Calcium.miss.modulationIdxM(Calcium.missNeuron),'k.'),hold on
%     for i = Calcium.missNeuron
%         line([1 2],[Calcium.hit.modulationIdxM(i) Calcium.miss.modulationIdxM(i)]),xlim([0.5 2.5])
%     end
% 

blah = mean(Calcium.hit.DeltaFoverF,3);
t = max(blah,[],2);
[~,idx] =  sort(t);
figure,imagesc(mean(Calcium.hit.motorNeuronTrace,3)),colormap(hot)
caxis([0.1 .6]),title('Motor Neuron')

figure,imagesc([mean(Calcium.hit.hitNeuronTrace,3);]),colormap(hot)
caxis([0.1 .6]),title('Hit Neuron')

figure,imagesc([mean(Calcium.hit.missNeuronTrace,3);mean(Calcium.miss.missNeuronTrace,3)]),colormap(hot)
caxis([0.1 .6]),title('Miss Neuron')
end
end