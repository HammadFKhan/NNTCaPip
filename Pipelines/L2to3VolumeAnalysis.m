

for ii = 1:floor(30168/603);
    idxAdjust = (ii-1)*603;
    stimWin = 100+idxAdjust:350+idxAdjust; %touch epoch frame array, TODO:verify via analog trace
    % Now we calculate the touch evoked event for all cells for a single trial
    for i = 1:size(DeltaFoverF,1) %num cells
        touchcalcium_avg{i,ii} = STA(DeltaFoverF(i,stimWin),2,100);%std, window (frames)
    end
    stimWin = 350+idxAdjust:603+idxAdjust;
    % now do it again but for the whole trial so we can compare
    for i = 1:size(DeltaFoverF,1) %num cells
        nontouchcalcium_avg{i,ii} = STA(DeltaFoverF(i,stimWin),2,100);%std, window (frames)
    end
end

% Generate a metric for statistical comparison
touchEvoked = cellfun(@(x) size(x,2), touchcalcium_avg); %number of touch evoked calcium event
allEvoked = cellfun(@(x) size(x,2), allcalcium_avg); % number of all calcium events
touchModulation = (touchEvoked)/(allEvoked); % lets assume the null would be 0.5