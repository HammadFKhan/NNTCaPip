function Calcium = leverCaModulation(DeltaFoverF,Spikes,Behaviour)
%Genenerate behavior triggered calcium response
caM = mean(DeltaFoverF,2);
for i = 1:Behaviour.nHit
    Calcium.hit.caTrace(i).DeltaFoverF = DeltaFoverF(:,Behaviour.hitTrace(i).caIndex);
    Calcium.hit.modulationIdx(:,i) = single(1.2*mean(Calcium.hit.caTrace(i).DeltaFoverF,2)>caM);
    Calcium.hit.modulationIdxM = mean(Calcium.hit.modulationIdx,2);
    Calcium.hit.caTraceCombined(:,:,i) = Calcium.hit.caTrace(i).DeltaFoverF;
    Calcium.hit.spikeCombined(:,:,i) = Spikes(:,Behaviour.hitTrace(i).caIndex);
end
for i = 1:Behaviour.nMiss
    Calcium.miss.caTrace(i).DeltaFoverF = DeltaFoverF(:,Behaviour.missTrace(i).caIndex);
    Calcium.miss.modulationIdx(:,i) = single(1.2*mean(Calcium.miss.caTrace(i).DeltaFoverF,2)>caM);
    Calcium.miss.modulationIdxM = mean(Calcium.miss.modulationIdx,2);
    Calcium.miss.caTraceCombined(:,:,i) = Calcium.hit.caTrace(i).DeltaFoverF;
    Calcium.miss.spikeCombined(:,:,i) = Spikes(:,Behaviour.missTrace(i).caIndex);
end

% Define motor neurons (modulation during either hit/miss trials)
Calcium.motorNeuron = unique(find(Calcium.hit.modulationIdxM>0.5 | Calcium.miss.modulationIdxM>0.5));
% Define hit/reward neurons (modulated stronger during hit trials)
Calcium.hitNeuron = unique(find(Calcium.hit.modulationIdxM>0.5 & Calcium.hit.modulationIdxM>Calcium.miss.modulationIdxM));
% Define miss neurons (modulated stronger during miss trials)
Calcium.missNeuron = unique(find(Calcium.miss.modulationIdxM>0.5 & Calcium.miss.modulationIdxM>Calcium.hit.modulationIdxM));

% Build lever task calcium traces
for i = 1:Behaviour.nHit
    % Do for motor neuron
    [~,caMaxIdx] = max(Calcium.hit.caTrace(i).DeltaFoverF(Calcium.motorNeuron,:),[],2);
    [~,idx] = sort(caMaxIdx);
    motorNeuronTrace = Calcium.hit.caTrace(i).DeltaFoverF(Calcium.motorNeuron,:);
    Calcium.hit.motorNeuronTrace(:,:,i) = motorNeuronTrace(idx,:);
    
    % Do for hit neuron
    [~,caMaxIdx] = max(Calcium.hit.caTrace(i).DeltaFoverF(Calcium.hitNeuron,:),[],2);
    [~,idx] = sort(caMaxIdx);
    hitNeuronTrace = Calcium.hit.caTrace(i).DeltaFoverF(Calcium.hitNeuron,:);
    Calcium.hit.hitNeuronTrace(:,:,i) = hitNeuronTrace(idx,:);
    
    % Do for miss neuron
    [~,caMaxIdx] = max(Calcium.hit.caTrace(i).DeltaFoverF(Calcium.missNeuron,:),[],2);
    [~,idx] = sort(caMaxIdx);
    missNeuronTrace = Calcium.hit.caTrace(i).DeltaFoverF(Calcium.missNeuron,:);
    Calcium.hit.missNeuronTrace(:,:,i) = missNeuronTrace(idx,:);
end
for i = 1:Behaviour.nMiss
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

%% plot some stuff
plotOn = 1;
if plotOn
    %figure,boxplot([Calcium.hit.modulationIdxM,Calcium.miss.modulationIdxM])
    figure,plot(ones(1,size(caM,1)),Calcium.hit.modulationIdxM,'k.'),hold on,title('All modulated neurons')
    plot(2*ones(1,size(caM,1)),Calcium.miss.modulationIdxM,'k.'),hold on
    for i = 1:size(caM,1)
        line([1 2],[Calcium.hit.modulationIdxM(i) Calcium.miss.modulationIdxM(i)]),xlim([0.5 2.5])
    end
    
    figure,plot(ones(1,size(Calcium.motorNeuron,1)),Calcium.hit.modulationIdxM(Calcium.motorNeuron),'k.')...
        ,hold on,title('Motor Neuron')
    plot(2*ones(1,size(Calcium.motorNeuron,1)),Calcium.miss.modulationIdxM(Calcium.motorNeuron),'k.'),hold on
    for i = Calcium.motorNeuron
        line([1 2],[Calcium.hit.modulationIdxM(i) Calcium.miss.modulationIdxM(i)]),xlim([0.5 2.5])
    end
    
    figure,plot(ones(1,size(Calcium.hitNeuron,1)),Calcium.hit.modulationIdxM(Calcium.hitNeuron),'k.')...
        ,hold on,title('Hit Responsive Neuron')
    plot(2*ones(1,size(Calcium.hitNeuron,1)),Calcium.miss.modulationIdxM(Calcium.hitNeuron),'k.'),hold on
    for i = Calcium.hitNeuron
        line([1 2],[Calcium.hit.modulationIdxM(i) Calcium.miss.modulationIdxM(i)]),xlim([0.5 2.5])
    end
    
    figure,plot(ones(1,size(Calcium.missNeuron,1)),Calcium.hit.modulationIdxM(Calcium.missNeuron),'k.')...
        ,hold on,title('Miss Responsive Neuron')
    plot(2*ones(1,size(Calcium.missNeuron,1)),Calcium.miss.modulationIdxM(Calcium.missNeuron),'k.'),hold on
    for i = Calcium.missNeuron
        line([1 2],[Calcium.hit.modulationIdxM(i) Calcium.miss.modulationIdxM(i)]),xlim([0.5 2.5])
    end 
    
    figure,imagesc([mean(Calcium.hit.motorNeuronTrace,3);mean(Calcium.miss.motorNeuronTrace,3)]),colormap(hot)
    caxis([0.1 .6]),title('Motor Neuron')
    
    figure,imagesc([mean(Calcium.hit.hitNeuronTrace,3);mean(Calcium.miss.hitNeuronTrace,3)]),colormap(hot)
    caxis([0.1 .6]),title('Hit Neuron')
    
    figure,imagesc([mean(Calcium.hit.missNeuronTrace,3);mean(Calcium.miss.missNeuronTrace,3)]),colormap(hot)
    caxis([0.1 .6]),title('Miss Neuron')
end
end