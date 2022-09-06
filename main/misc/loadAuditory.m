function [y,Fs] = loadAuditory(y,Fs)
if nargin==0
    
    % Load auditory files and do some baseline analysis
    [filename, pathname] = uigetfile({'*.wav;*.wav'}, 'Pick a image video file');
    if isequal(filename,0) || isequal(pathname,0)
        disp('User pressed cancel')
        return
    else
        disp(['User selected ', fullfile(pathname, filename)])
    end
    [y,Fs] = audioread(fullfile(pathname, filename));
    
else
    % fft
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
    % PSD
    disp('Plotting audio PSD...')
    [Pxx,w] = periodogram(y,[],[24000],Fs);
    figure,plot(w,10*log10(Pxx)), hold on
    xlim([500 15000])
    xlabel('Frequency'),ylabel('Power dB');
    
    
    % Wavelet
    disp('Plotting audio wavelet...')
    [wavelet,f] = cwt(y,Fs,'FrequencyLimit',[500 15000]);
    figure,imagesc(0:length(y)/Fs,f,abs(wavelet)),axis xy,colormap(jet);
end