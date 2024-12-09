%% Ensemble Stability Anova stats (SI included)
load('F:\DataD_Drive_backup\Ensemble_redo\All\CTX_stats.mat')
load('F:\DataD_Drive_backup\Ensemble_redo\All\STR_stats.mat')
dat = totalMaxRecruitment;
temp = dat(:);
temp(isnan(temp)) = [];
group = "";
mouse = "";
mouseNum = 3;
for n = 1:6
    group = [group;repmat(num2str(n-1),length(dat(~isnan(dat(:,n)),1)),1)];
    mLen = floor(length(dat(~isnan(dat(:,n)),1))/mouseNum);
    for nMouse = 1:mouseNum-1
        mouse = [mouse;repmat(num2str(nMouse),mLen,1)];
    end
    blah = length(dat(~isnan(dat(:,n)),1))-(mouseNum-1)*mLen;
    mouse = [mouse;repmat(num2str(nMouse+1),blah,1)];
end
group(1) = [];
mouse(1) = [];
[p,t,stats] = anovan(temp,{group mouse},'model','interaction','varnames',{'group','mouse'});
%% Pairwise stats
load('Y:\Hammad\Ephys\PFFProject\Stats\CTXPairwiseStats')
dat = CTXpair;
temp = dat(:);
temp(isnan(temp)) = [];
groupPair = "";
mousePair = "";
mouseNum = 3;
for n = 1:6
    groupPair = [groupPair;repmat(num2str(n-1),length(dat(~isnan(dat(:,n)),1)),1)];
    mLen = floor(length(dat(~isnan(dat(:,n)),1))/mouseNum);
    for nMouse = 1:mouseNum-1
        mousePair = [mousePair;repmat(num2str(nMouse),mLen,1)];
    end
    blah = length(dat(~isnan(dat(:,n)),1))-(mouseNum-1)*mLen;
    mousePair = [mousePair;repmat(num2str(nMouse+1),blah,1)];
end
groupPair(1) = [];
mousePair(1) = [];
[p,t,stats] = anovan(temp,{groupPair mousePair},'model','interaction','varnames',{'group','mouse'});
multcompare(stats);
%% Ensemble Num
load('F:\DataD_Drive_backup\Ensemble_redo\All\CTX_stats.mat')
%load('F:\DataD_Drive_backup\Ensemble_redo\All\STR_stats.mat')
dat = totalEnsembleSize;
temp = dat(:);
temp(isnan(temp)) = [];
group = "";
mouse = "";
mouseNum = 3;
for n = 1:6
    group = [group;repmat(num2str(n-1),length(dat(~isnan(dat(:,n)),1)),1)];
    mLen = floor(length(dat(~isnan(dat(:,n)),1))/mouseNum);
    for nMouse = 1:mouseNum-1
        mouse = [mouse;repmat(num2str(nMouse),mLen,1)];
    end
    blah = length(dat(~isnan(dat(:,n)),1))-(mouseNum-1)*mLen;
    mouse = [mouse;repmat(num2str(nMouse+1),blah,1)];
end
group(1) = [];
mouse(1) = [];
ensemble = (1:length(temp))';
nest_matrix = zeros(3,3); % matrix that specify which factors are nested
nest_matrix(2,1) = 1; % mouse is nested in the treatment groups
nest_matrix(3,2) = 1; % neuron is nested in mouse
[p1,tl,stats1,terms1] = anovan(temp, {group, mouse, ensemble}, 'varnames', {'treatment group', 'mouse', 'ensemble'}, 'nested', nest_matrix, 'random', [2,3]);
