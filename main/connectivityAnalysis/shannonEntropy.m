function informationEntropy = shannonEntropy(rankEnsembles,populationSize)
%% Shannon Entropy Analysis
% Probability of each ensemble
ensembleSize = cellfun(@size,rankEnsembles,'UniformOutput',false)';
ensembleSize = cell2mat(ensembleSize);
ensembleSize(:,2) = [];
totSize = sum(ensembleSize);
% P_c = zeros(length(ensembleSize),1);
% P_e = zeros(length(ensembleSize),1);
% I_c = zeros(length(ensembleSize),1);

P_c = ensembleSize/totSize; %Probability of ensemble based on cell ensemble size
P_s = ensembleSize/populationSize; % Probability based on the actual size of the field of view
P_e = repmat(length(ensembleSize)^-1,length(ensembleSize),1); %Probability of ensemble activation based on # detected
I_c = -log2(P_c); %Ensemble self information
% Extend to general probability distribution for information entropy of
% recording session
H_c = sum(-P_c.*log2(P_c));
disp(['Information Entropy of Ensemble Activation: ' num2str(H_c) ' bits']); %Distribution of cells making up ensembles
H_e = sum(-P_e.*log2(P_e)); %Self information and entropy are equal since the probabilities are the same. 
disp(['Information Entropy of Ensemble Frequency: ' num2str(H_e) ' bits']); %Distrubution of ensemble numbers
H_s = sum(-P_s.*log2(P_s)); %Self information and entropy are equal since the probabilities are the same. 
disp(['Information Entropy of Ensemble Size: ' num2str(H_s) ' bits']); %Distrubution of ensemble numbers
% Output
informationEntropy.P_c = P_c;
informationEntropy.P_s = P_s;
informationEntropy.P_e = P_e;
informationEntropy.I_c = I_c;
informationEntropy.H_c = H_c;
informationEntropy.H_e = H_e;
informationEntropy.H_s = H_s;

end