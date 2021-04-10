%% EGB242 Lab - Quantization 3
clear all, close all, clc

load('EGB242_LabC.mat')

% Choose music - 16 bit 44.1k acting as our Continuous-Time Continuous-Amplitude signal
sig = s2;
% sound(sig,fs) % Comment out when done
t = linspace(0,length(sig)/fs,length(sig)+1); t = t(1:end-1);
figure(1), plot(t,sig), axis('tight'), hold on

% Uniform quantizer - assume normalised limits (+/-1)
maxLim = 1; % max(sig)
minLim = -1; % min(sig)
bits = 8;  numLevs = 2^bits; % Start at 8 then 6,4,3
delta = (maxLim - minLim)/numLevs;

% Sample every I samples - Discrete-Time Continuous-Amplitude
I = 10;
fsSamp = fs/I;
sigSamp = sig(1:I:end);
tSamp = linspace(0,length(sigSamp)/fsSamp,length(sigSamp)+1);
tSamp = tSamp(1:end-1);
plot(tSamp,sigSamp,'r*',tSamp,sigSamp,'m--')
% axis([4.274 4.284 -1 1]) % manually zoom in to part with poorly sampled high frequency content

% sound(sigSamp,fsSamp) % Comment out when done
% Note - fewer high frequencies, some distortion/errors coming from being unable
% to replicate frequencies above Nyquist rate
% --> Aliasing (very obvious in s4) - pre-filter to avoid this (LPF with cutoff below Nyquist
% frequency

% Quantize (Uniform MR) - Discrete-Time Discrete-Amplitude - Digital
sigQuant = delta*(floor(sigSamp/delta) + 1/2);

plot(tSamp,sigQuant,'go')

% sound(sigQuant,fsSamp) % Comment out when done
% Note - Quantization noise can be heard as hiss in the background
% Try smaller bit depths to increase effect

% Frequency Response comparison
dF = fs/length(t);
f = -fs/2:dF:fs/2-dF;

fSamp = linspace(-fsSamp/2,fsSamp/2,length(tSamp)+1);
fSamp = fSamp(1:end-1);

figure(2), subplot(3,1,1), plot(f,abs(fftshift(fft(sig)/fs)),'b')
axis([fSamp(1) fSamp(end) 0 0.2])
subplot(3,1,2), plot(fSamp,abs(fftshift(fft(sigSamp)/fsSamp)),'m')
subplot(3,1,3), plot(fSamp,abs(fftshift(fft(sigQuant)/fsSamp)),'g')