% set general variables
Fs = 44100;  % sample frequency
nf = Fs / 2; % nyquist frequency
d = 2.0;     % duration (time)
n = Fs * d;  % number of samples
nh = n / 2;  % half number of samples
N = 3;      % Harmonic Degree
F0 = 2000;
%% Generate Harmonics
secs = 2;
t  = linspace(0, secs, Fs*secs+1);              % Time Vector + 1 sample
t(end) = [];  % remove extra sample
harmonic = F0:F0:F0*N;
st = zeros(1,length(t));
for i = 1:N
w = 2*pi*harmonic(i);
s(i,:) = sin(w*t);
end

st = conv(s(1,:),s(2,:),'same');
st = conv(st,s(3,:),'same');
%%
for i = 1:3
    sound([st], Fs)                              % Produce Tone then second Tone As Sound
    disp(['Sound sequence: ' num2str(i)])
    pause(secs*2)
end
%%
y = st;
disp('Plotting audio FFT')
Y = fft(y);
L = length(y);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
figure,
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

disp('Plotting audio PSD...')
[Pxx,w] = periodogram(y,[],[24000],Fs);
figure,plot(w,10*log10(Pxx)), hold on
xlim([500 15000])
xlabel('Frequency'),ylabel('Power dB');