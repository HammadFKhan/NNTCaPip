function batchSpikes = TrialByTrial(batchData)

% ### Requires Batch Structure ###
cellSpikes = struct2cell(batchData(:));
spikeSize = zeros(1,size(batchData,2));
for i = 1:size(batchData,2)
    spikeSize(i) = size(cellSpikes{1,i},1);
end

n = size(batchData,2); % number of matrices
r = min(spikeSize); % row of matrices
c = size(cellSpikes{1,i},2); % column of matrices
batchSpikes = zeros(r,c,n);
for ii = 1:n
    spikes = cellSpikes{1,ii};
    try
    batchSpikes(:,:,ii) = spikes(1:r,:);
    catch ME
        warning('Time dimensions do not agree.');
        disp('Not all trials were analyzed!')
        break
    end
end

batchSpikes = reshape(batchSpikes,r,c*n); % Reshape into single matrix
end