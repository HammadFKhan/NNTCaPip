load('Y:\Hammad\Ephys\PFFProject\Stats\SpikeRateCTXAnova.mat')
%%
temp = spikeDat(:,1:6);
temp = temp(:);
temp(isnan(temp)) = [];
group = "";
mouse = "";
mouseNum = 4;
for n = 1:6
    group = [group;repmat(num2str(n-1),length(spikeDat(~isnan(spikeDat(:,n)),1)),1)];
    mLen = floor(length(spikeDat(~isnan(spikeDat(:,n)),1))/mouseNum);
    for nMouse = 1:mouseNum-1
        mouse = [mouse;repmat(num2str(nMouse),mLen,1)];
    end
    mouse = [mouse;repmat(num2str(nMouse+1),length(spikeDat(~isnan(spikeDat(:,n)),1))-(mouseNum-1)*mLen,1)];
end
group(1,:) = [];
mouse(1,:) = [];
[p,t,stats] = anovan(temp,{group mouse},'model','interaction','varnames',{'group','mouse'});
[results,~,~,gnames] = multcompare(stats,"Dimension",[1 1]);
%%
tbl1 = array2table(results,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])
%% 
figure,customBoxplot(spikeDat)
%%
neuronMouse = [];
for n= 1:6
    idx = strcmp(group,num2str(n-1));
    for nn = 1:4
        neuronMouse(nn,n) = mean(temp(idx & strcmp(mouse,num2str(nn))));
    end
end
%%
[p,~,stats] = anova2(neuronMouse,2)
c = multcompare(stats,'estimate','column')
