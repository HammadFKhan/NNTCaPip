%% Proteomic analysis

cg = clustergram([CTXMon CTXPFF(:,[1,4:6])],'Rowlabels',clabels,'standardize','row');
cg.Colormap = redbluecmap

cg = clustergram([STRMon STRPFF],'Rowlabels',clabels,'standardize','row');
cg.Colormap = redbluecmap
%%
figure,scatter(strRegulated(:,1),strRegulated(:,2),10,'filled'), hold on
scatter(strUpRegulated(:,1),strUpRegulated(:,2),10,'filled')
scatter(strDownRegulated(:,1),strDownRegulated(:,2),10,'filled')
%% Pie Plot for Q-value hits
ctxPie = [size(ctxUpRegulated,1), size(ctxDownRegulated,1),size(ctxRegulated,1)];
figure,pie(ctxPie),title('CTX')
strPie = [size(strUpRegulated,1), size(strDownRegulated,1),size(strRegulated,1)];
figure,pie(strPie),title('STR')
%% NMDA
ctxNMDAlabels = categorical(ctxNMDAlabel);
figure,bar(ctxNMDAlabels, ctxNMDA)
box off
set(gca,'TickDir','out')
strNMDAlabels = categorical(strNMDAlabel);
figure,bar(strNMDAlabels, strNMDA)
box off
set(gca,'TickDir','out')
%%
CTXMonNMDA = NMDA(:,1:5);
STRMonNMDA = NMDA(:,6:10);
CTXPFFNMDA = NMDA(:,11:16);
STRPFFNMDA = NMDA(:,17:22);
cg = clustergram([CTXMonNMDA(:,[1 2 3 4 5]) CTXPFFNMDA(:,2:6)],'Rowlabels',NMDAlabels,'standardize','row');
cg.Colormap = redbluecmap

cg = clustergram([STRMonNMDA STRPFFNMDA(:,1:5)],'Rowlabels',NMDAlabels,'standardize','row');
cg.Colormap = redbluecmap