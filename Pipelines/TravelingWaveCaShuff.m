%% Shuffled trials and check similarity
lateSpikeTrialsShuf = [];
nolateSpikeTrialsShuf = [];
surrogate = 2;
for i = 1:surrogate
    lateSpikeTrialsShuftemp = randi([1 45]...
        ,1,length(trialData.responsiveTrials.lateSpikeTrials));
    
    nolateSpikeTrialsShuftemp = randi([1 45]...
        ,1,length(trialData.responsiveTrials.noLateSpikeTrials));
    
    lateSpikeTrialsShuf = horzcat(lateSpikeTrialsShuf,lateSpikeTrialsShuftemp);
    nolateSpikeTrialsShuf = horzcat(nolateSpikeTrialsShuf,nolateSpikeTrialsShuftemp);
end
%%
% simM = [];
% Connected_ROI = [];
% Ls = [];
% % figure,
% count = 1;
% for i = nolateSpikeTrialsShuf  %t1
%     corr = correlation_dice(spikeTrials{i});
%     Connected_ROI{count} = Connectivity_dice(corr,ROI,0.4);
%     [NumActiveNodes{count},NodeList{count},NumNodes{count},NumEdges{count},SpatialCentroid{count},SpatialCentroidVariance{count},ActivityCentroid{count},ActivityCentroidVariance{count}, ActivityCoords{count}]...
%         = Network_Analysis(ROIcentroid,Connected_ROI{count});
%     [vectorized,sim_index] = cosine_similarity(spikeTrials{i},10);
% %     subplot(5,4,count),imagesc(sim_index),colormap(jet)
%     Ls(count) = mean(sim_index(sim_index>0.1),'all'); %mean center
%     count = count+1;
% end
% % t1 = vertcat(simValue{:});
%  figure,boxplot(Ls);ylim([0.0 0.5])
 %%
 %% Quantify Node Reactivation
lateSpikeCritNodes = findCritNode(lateSpikeTrialsShuf,ROI,spikeTrials,Spikes);
nolateSpikeCritNodes = findCritNode(nolateSpikeTrialsShuf,ROI,spikeTrials,Spikes);

%%
count = 1;
nls_ls = []
Connected_ROI = [];
for i = lateSpikeTrialsShuf%t2
    corr = correlation_dice(spikeTrials{i}(nolateSpikeCritNodes,:));
    Connected_ROI{count} = Connectivity_dice(corr,ROI,0.3);
    count = count+1;
end
nls_ls = vertcat(Connected_ROI{:});
figure,
Cell_Map_Dice(AverageImage,nls_ls,ROIcentroid(nolateSpikeCritNodes,:),4,1)
%%
l1 = [length(ls_ls) length(ls_nls) length(nls_nls) length(nls_ls)];
% simS = [ls_ls(:,3) ls_nls(:,3) nls_nls(:,3) nls_ls(:,3)];
figure,boxplot([l(1) l1(1)]);
%% Do a bunch
lateSpikeTrialsShuf = [];
nolateSpikeTrialsShuf = [];
surrogate = 10;
for ii = 1:surrogate
    disp(['Iteration: ' num2str(ii)]);
    lateSpikeTrialsShuf = randi([1 45]...
        ,1,length(trialData.responsiveTrials.lateSpikeTrials));
    
    nolateSpikeTrialsShuf = randi([1 45]...
        ,1,length(trialData.responsiveTrials.noLateSpikeTrials));
    
    lateSpikeCritNodes = findCritNode(lateSpikeTrialsShuf,ROI,spikeTrials,Spikes);
    nolateSpikeCritNodes = findCritNode(nolateSpikeTrialsShuf,ROI,spikeTrials,Spikes);
    
    % Ls:Ls
    ls_ls = []
    Connected_ROI = [];
    for i = lateSpikeTrialsShuf%t2
        corr = correlation_dice(spikeTrials{i}(lateSpikeCritNodes,:));
        Connected_ROI{count} = Connectivity_dice(corr,ROI,0.3);
        count = count+1;
    end
    ls_ls = vertcat(Connected_ROI{:});
    
    % Ls:nLs
    ls_nls = []
    Connected_ROI = [];
    for i = nolateSpikeTrialsShuf%t2
        corr = correlation_dice(spikeTrials{i}(lateSpikeCritNodes,:));
        Connected_ROI{count} = Connectivity_dice(corr,ROI,0.3);
        count = count+1;
    end
    ls_nls = vertcat(Connected_ROI{:});
    
    % nLs:nLs
    nls_nls = []
    Connected_ROI = [];
    for i = nolateSpikeTrialsShuf%t2
        corr = correlation_dice(spikeTrials{i}(nolateSpikeCritNodes,:));
        Connected_ROI{count} = Connectivity_dice(corr,ROI,0.3);
        count = count+1;
    end
    nls_nls = vertcat(Connected_ROI{:});
    
    % nLs:Ls
    nls_ls = []
    Connected_ROI = [];
    for i = lateSpikeTrialsShuf%t2
        corr = correlation_dice(spikeTrials{i}(nolateSpikeCritNodes,:));
        Connected_ROI{count} = Connectivity_dice(corr,ROI,0.3);
        count = count+1;
    end
    nls_ls = vertcat(Connected_ROI{:});
    l(:,ii) = [length(ls_ls) length(ls_nls) length(nls_nls) length(nls_ls)];
    simL(:,ii) = [mean(ls_ls(:,3)) mean(ls_nls(:,3)) mean(nls_nls(:,3)) mean(nls_ls(:,3))];
end
%%
figure,boxplot(simL(4,:)), box off, ylim([0 1])