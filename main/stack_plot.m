function stack_plot(DeltaFoverF,spacing,amplitude,normFlag)
if nargin < 4 || strcmp(normFlag,'')
    normFlag = 0;
end
if nargin < 2 || strcmp(amplitude,'')
    amplitude = 1;
end
if nargin < 2 || strcmp(spacing,'')
    spacing = 1;
end
% Fs = 30.048;
Fs = 10.5759;
x = length(DeltaFoverF(1,:));
y = length(DeltaFoverF(:,1));
DeltaFoverF = flip(DeltaFoverF);
time = 1/Fs:1/Fs:x/Fs;
[grad,~]=colorGradient([7 49 97]/255,[110 192 235]/255,y);
if normFlag
    for i = 1:y
        DeltaFoverF(i,:) = (DeltaFoverF(i,:)-min(DeltaFoverF(i,:)))/(max(DeltaFoverF(i,:))-min(DeltaFoverF(i,:)));
    end
end
baseline = spacing*max(DeltaFoverF,[],'all');
for i = 1:y
    gradient = i/y;
    plot(time,amplitude*DeltaFoverF(i,:)+(spacing*baseline),'LineWidth',1,'Color',grad(i,:)); hold on; % ,'Color',[.30 .835 .384]
    baseline = baseline + spacing*max(DeltaFoverF,[],'all');
    axis tight, box off,...
            set(gca,'YTick','','YTickLabel','');
end

end
