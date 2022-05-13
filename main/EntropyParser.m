
function [P_c,I_c,H_c,H_e,sizeE,sizeEdge] = EntropyParser(Entropy)
P_c = Entropy.informationEntropy.P_c;
I_c = Entropy.informationEntropy.I_c;
H_c = Entropy.informationEntropy.H_c;
H_e = Entropy.informationEntropy.H_e;

sizeE = cellfun(@size,Entropy.rankedEnsembles,'UniformOutput',false);
sizeE = cell2mat(sizeE);
sizeE(sizeE==1) = [];
sizeE = sizeE/max(sizeE);

sizeEdge = cell2mat(Entropy.rankedEdges);
sizeEdge = sizeEdge/max(sizeEdge);
