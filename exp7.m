clc;clear all;close all;
fs = 1e6;                   
numSamples = 10000;          
numPaths = 5;                
maxDelay = 3e-6;             
dopplerShift = 100;          
impulseSignal = [1; zeros(numSamples-1, 1)]; 
rayleighChan = comm.RayleighChannel( ...
    'SampleRate', fs, ...
    'PathDelays', linspace(0, maxDelay, numPaths), ... 
    'AveragePathGains', [-2 -3 -6 -8 -10], ...          
    'MaximumDopplerShift', dopplerShift, ...            
    'NormalizePathGains', true);
rxImpulseSignal = rayleighChan(impulseSignal);
timeAxis = (0:numSamples-1)/fs; 
figure;
stem(timeAxis(1:100), 20*log10(abs(rxImpulseSignal(1:100))));  
title('Impulse Response of Frequency-Selective Rayleigh Fading Channel');
xlabel('Time (s)');ylabel('Gain (dB)');grid on;
NFFT = 1024;  
freqResponse = fft(rxImpulseSignal, NFFT); 
freqAxis = linspace(-fs/2, fs/2, NFFT);  
figure;
plot(freqAxis/1e6, 20*log10(abs(fftshift(freqResponse)))); 
title('Frequency Response of Frequency-Selective Rayleigh Fading Channel');
xlabel('Frequency (MHz)');ylabel('Magnitude (dB)');grid on;

