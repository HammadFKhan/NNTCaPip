function LFP = Ca_LFP(Catime,loadFlag,out)
plotFig = 0;
if loadFlag
    out=import_wcp;
end
Vm = out.S{1};
Im = out.S{2};
T = out.T;
dt = out.T(2)-out.T(1);% time
Fs = 1/dt; 

% Note that under the digidata 1440 using WinWCP the
% sampling rate of the software does not seem to match
% the sampling rate of the actual data. To check this we 
% can compare using the behavior/calcium imaging data since
% time is constant

% output
LFP.out = out;

%% Vm
% Fix DC offset
VmTotal = [];
for i = 1:size(Vm,2)
    VmTotal = [VmTotal;Vm(:,i)];
end

VmTrigAlign = behaviorTrig(out,VmTotal);
% Adjust units to uV for detection analysis
VmTrigAlign = VmTrigAlign*1000;

% Compare Fs rate
CaFs = length(VmTrigAlign)/Catime(end);
if CaFs ~= Fs
    Fs = CaFs;
    warning('Digidata Fs has been matched to calculated calcium Fs!')
end

%% Im
% Fix DC offset
ImTotal = [];
Im = Im-min(Im,[],'all');
for i = 1:size(Im,2)
    ImTotal = [ImTotal;Im(:,i)];
end
%% Broadband filter
Fc = 3;
Wn = Fc./(Fs/2);
b1 = fir1(10000,Wn,'high');
disp('Filtering...')
tic;d = filtfilt(b1,1,VmTrigAlign);toc;

Fc = [4 250];
Wn = Fc./(Fs/2);
b = fir1(10000,Wn,'bandpass');
tic;Vmfilt = filtfilt(b,1,d);toc;
%%
subThreshold = Vmfilt;
thresh = find(Vmfilt>1000);
tri = diff(thresh);
trig = find(tri~=1);
COM = thresh(trig+1,1);

for i = 1:size(COM,2)
    thresh = find(Vmfilt>1000);
    tri = diff(thresh);
    trig = find(tri~=1);
    COM = thresh(trig+1,1);
    subThreshold(COM(1)-(0.01*Fs):COM(1)+(0.02*Fs)) = 500;
end
%% Filter
filtered_data = customFilt(Vmfilt',Fs,[10 30]);
Beta = IntrabetaBurstDetection(filtered_data',Fs);
%% Output
if plotFig
    figure,plot(0:dt:(length(VmTotal)-1)/Fs,VmTotal);xlim([0 length(VmTotal)/Fs])
    figure,plot(0:dt:(length(ImTotal)-1)/Fs,ImTotal);xlim([0 length(ImTotal)/Fs])
end
LFP.beta = Beta;
LFP.Vmfilt = Vmfilt;
LFP.betaLFP = filtered_data;
LFP.Fs = Fs;
LFP.time = T;
LFP.dt = dt;

