%pathname = 'F:\DataD_Drive_backup\Ensemble_redo\CTX_Control';
% pathname = 'F:\DataD_Drive_backup\Ensemble_redo\STR_Control';
%pathname = 'F:\DataD_Drive_backup\Ensemble_redo\UninjectedGcamp';
filesnames = dir(fullfile(pathname,'*.mat'));
for n = 1:length(filesnames)
    load(fullfile(filesnames(n).folder,filesnames(n).name))
    Mice(n).Ensemble = Ensemble;
    temp = correlation_dice(Ensemble.ensemble);
    diceCorr{n} = temp(temp>0);
    sim_index{n} =(Ensemble.sim_index(Ensemble.sim_index>0));
end
%% Representative correlation map
cmap = cmocean('thermal');
id = 3;
correctionFactor = floor(size(Mice(id).Ensemble.ensemble,2)/5)*5;
[~,simIndex] = cosine_similarity(Mice(id).Ensemble.ensemble(:,1:correctionFactor),5);
figure,imagesc(simIndex);colormap(cmap);axis([0 200 0 200]),caxis([0 1]),colorbar
% dat = correlation_dice(Mice(id).Ensemble.ensemble);
% figure,imagesc(dat);colormap(hot);caxis([0 0.5])
%%
dcM = cellfun(@mean,diceCorr);
dcSE = cellfun(@(x) std(x)/sqrt(size(x,1)/10), diceCorr);
figure,
errorbar(1:length(dcM),dcM,dcSE);
xlim([0.5 length(dcM)+.5]),ylim([0 0.4])
box off
set(gca,'tickdir','out','fontsize',16')
ylabel('Pairwise Synchrony')
figure,customBoxplot(dcM')
box off
set(gca,'tickdir','out','fontsize',16')
ylabel('Pairwise Synchrony')
ylim([0 0.4])
%%
simM = cellfun(@mean,sim_index);
dcSE = cellfun(@(x) std(x)/sqrt(size(x,1)/1000), sim_index);
figure,
errorbar(1:length(simM),simM,dcSE);
xlim([0.5 5.5]),ylim([0 0.6])
box off
set(gca,'tickdir','out','fontsize',16')
ylabel('Ensemble Stability')
%%
figure,customBoxplot(CTXMonStab)
box off
set(gca,'tickdir','out','fontsize',16')
ylabel('Ensemble Stability')
ylim([0 0.4])
%%
load('F:\DataD_Drive_backup\Ensemble_redo\CTX_Control\totalData.mat')
% Combine pairwise plots for analysis
dat1 = CTXMonPair;
dat2 = STRMonPair;
dat3 = WTPairwise;
pairwisetot = zeros(max([length(dat1),length(dat2),length(dat3)]),3);
pairwisetot(1:length(dat1),1) = dat1;
pairwisetot(1:length(dat2),2) = dat2;
pairwisetot(1:length(dat3),3) = dat3;
pairwisetot(pairwisetot==0) = NaN;

dat1 = CTXMonStab;
dat2 = STRMonStab;
dat3 = WTStability;
ensemblestabtot = zeros(max([length(dat1),length(dat2),length(dat3)]),3);
ensemblestabtot(1:length(dat1),1) = dat1;
ensemblestabtot(1:length(dat2),2) = dat2;
ensemblestabtot(1:length(dat3),3) = dat3;
ensemblestabtot(ensemblestabtot==0) = NaN;

anova1(pairwisetot)
anova1(ensemblestabtot);
%% W2W12 FOV
load("F:\DataD_Drive_backup\ControlData_ExtendedFig\RecordM1_Big_00003-2")
Spikes = rasterizeDFoF(dDeltaFoverF,2.5,0.01);
FOV1 = correlation_dice(Spikes);
load("F:\DataD_Drive_backup\ControlData_ExtendedFig\M1_ROI1_00001-1")
FOV2 = correlation_dice(Spikes);
figure,histogram(FOV1(FOV1>0),0:0.01:1,'normalization','probability','edgecolor','none'),hold on
histogram(FOV2(FOV2>0),0:0.01:1,'normalization','probability','edgecolor','none'),hold on
xlim([0 0.2])