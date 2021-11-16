%% Assignment 2, Part 2 (EEG signal analysis)
%  Do not change before line 35
%  You will need to have generated A2P2Data.mat from 
%  GenerateAssignment2Data.m before working with this file.

%  Clearing and preparing the workspace
clear; clc; close all;

%  Load assignment data from A2P1Data.mat.
load('A2P2Data.mat');  

%=================================================================%
%
% Refer to the assignment sheet for details.
% Names of the variables are important,
% e.g. 'a1' is considered a different variable to 'A1'.
% Make sure variables have been declared as they appear in the brief.
%
% Ts - sampling period
% t - time domain vector
% MUXSIG - Fourier Transform
% k - frequency vector
% fshift - frequency shifts
% Mag - magnitude
% Phase- phase
% xdm - EEG data
% XDM - Fourier Transform
% freqResponse - frequency response of systems
% impresp - impulse response of chosen system
% EEG - filtered signals
% eeg - time domain of filtered signals
% Conveq - convolution
%
%====Enter your code below this line================================


%% Question 2.1 %% Represent the spectrum analyzer 

%determing sampling period%
Ts = 1/fs;

%creating time vector%
t = linspace(0,Ts*length(muxSignal),length(muxSignal)+1);
t(end) = [];

%plotting muxSignal vs. time%
figure(1)
plot(t, muxSignal);
xlabel('Time [s]');
ylabel('Amplitude');
title('muxSignal time domain spectrum');

%computing Fourier transform of muxSignal%

MUXSIG = fft(muxSignal);

%creating frequency vector k%
k = linspace(-fs/2,fs/2, length(MUXSIG)+1);
k(end) = [];

%plotting the magnitude spectra of MUXSIG%

figure(2)
plot(k, abs(MUXSIG));
xlabel('Frequency [Hz]');
ylabel('Magnitude');
title('MUXSIG frequency spectrum');

%% Question 2.2%% determining demultiplexing parameters

%Shifting and scaling MUXSIG%
MUXSIG = fftshift(MUXSIG)/fs;

%identifying frequency shifts%
%noting peak heights above 0.06 magnitude%

[PKS,LOCS] = findpeaks(abs(MUXSIG), 'MinPeakHeight', 0.06);
LOCS = LOCS(6:end);

%storing frequency shifts in vector 'fshift'%
fshift = k(LOCS);

%finding corresponding magnitude and phase%
Mag = abs(MUXSIG(LOCS))';
Phase = angle(MUXSIG(LOCS))';

%transpose muxSignal to row vector%
muxSignal_tpose = muxSignal';

%removing freq shifts with xdm function%
[xdm] = FDMDemux(muxSignal_tpose,t,Mag,fshift,Phase);

%% Question 2.4%% reviewing individual signals%

XDM = zeros(5, length(xdm));

XDM(1,:) = fft(xdm(1,:));
XDM(2,:) = fft(xdm(2,:));
XDM(3,:) = fft(xdm(3,:));
XDM(4,:) = fft(xdm(4,:));
XDM(5,:) = fft(xdm(5,:));

figure(3)
subplot(3,2,1)
plot(k, abs(XDM(1,:)))
xlabel('Frequency [Hz]')
ylabel('Magnitude')
title('EEG signal astronaut 1')
subplot(3,2,2)
plot(k, abs(XDM(2,:)))
xlabel('Frequency [Hz]')
ylabel('Magnitude')
title('EEG signal astronaut 2')
subplot(3,2,3)
plot(k, abs(XDM(3,:)))
xlabel('Frequency [Hz]')
ylabel('Magnitude')
title('EEG signal astronaut 3')
subplot(3,2,4)
plot(k, abs(XDM(4,:)))
xlabel('Frequency [Hz]')
ylabel('Magnitude')
title('EEG signal astronaut 4')
subplot(3,2,5)
plot(k, abs(XDM(5,:)))
xlabel('Frequency [Hz]')
ylabel('Magnitude')
title('EEG signal astronaut 5')

XDM_shift = zeros(5, length(XDM));

XDM_shift(1,:) = fftshift(XDM(1,:))/fs;
XDM_shift(2,:) = fftshift(XDM(2,:))/fs;
XDM_shift(3,:) = fftshift(XDM(3,:))/fs;
XDM_shift(4,:) = fftshift(XDM(4,:))/fs;
XDM_shift(5,:) = fftshift(XDM(5,:))/fs;

figure(4)
subplot(3,2,1)
plot(k, abs(XDM_shift(1,:)))
xlabel('Frequency [Hz]')
ylabel('Magnitude')
title('EEG signal astronaut 1')
subplot(3,2,2)
plot(k, abs(XDM_shift(2,:)))
xlabel('Frequency [Hz]')
ylabel('Magnitude')
title('EEG signal astronaut 2')
subplot(3,2,3)
plot(k, abs(XDM_shift(3,:)))
xlabel('Frequency [Hz]')
ylabel('Magnitude')
title('EEG signal astronaut 3')
subplot(3,2,4)
plot(k, abs(XDM_shift(4,:)))
xlabel('Frequency [Hz]')
ylabel('Magnitude')
title('EEG signal astronaut 4')
subplot(3,2,5)
plot(k, abs(XDM_shift(5,:)))
xlabel('Frequency [Hz]')
ylabel('Magnitude')
title('EEG signal astronaut 5')

%% Question 2.5, 2.6 & 2.7 evaluating non ideal filters

factorTF(sys(1))

factorTF(sys(2))

factorTF(sys(3))

factorTF(sys(4))

ltiview(sys(3))
%all filters have poles on negative side of cartesian plane "stable"
%filter 1 takes a long time to settle, bode plot off also
%filter 2 step response doesnt stabilize at 1, low pass filter
%filter 3 settles well, high pass judging by bode plot and is stable
%filter 4 doesnt stablize at 1 for step response, low pass filter 

%% Question 2.8 %% applying non ideal filter

factorTF(sys(3))

%impulse response of system
h_f = impulse(sys(3),t)';

H_f = fft(h_f);
H_f = fftshift(H_f)/fs;

%multiplication of impulse resp 
EEG(1,:) = H_f.*XDM_shift(1,:);
EEG(2,:) = H_f.*XDM_shift(2,:);
EEG(3,:) = H_f.*XDM_shift(3,:);
EEG(4,:) = H_f.*XDM_shift(4,:);
EEG(5,:) = H_f.*XDM_shift(5,:);

%plotting non-ideal filtered signals in freq domain
subplot(3,2,1)
plot(k, abs(EEG(1,:)))
xlabel('Frequency [Hz]')
ylabel('Magnitude')
title('Filtered EEG signal 1 - Frequency spectrum')
subplot(3,2,2)
plot(k, abs(EEG(2,:)))
xlabel('Frequency [Hz]')
ylabel('Magnitude')
title('Filtered EEG signal 2 - Frequency spectrum')
subplot(3,2,3)
plot(k, abs(EEG(3,:)))
xlabel('Frequency [Hz]')
ylabel('Magnitude')
title('Filtered EEG signal 3 - Frequency spectrum')
subplot(3,2,4)
plot(k, abs(EEG(4,:)))
xlabel('Frequency [Hz]')
ylabel('Magnitude')
title('Filtered EEG signal 4 - Frequency spectrum')
subplot(3,2,5)
plot(k, abs(EEG(5,:)))
xlabel('Frequency [Hz]')
ylabel('Magnitude')
title('Filtered EEG signal 5 - Frequency spectrum')



%bringing EEG to time domain
eeg(1,:) = ifftshift(EEG(1,:))*fs;
eeg(1,:) = ifft(eeg(1,:));
eeg(2,:) = ifftshift(EEG(2,:))*fs;
eeg(2,:) = ifft(eeg(2,:));
eeg(3,:) = ifftshift(EEG(3,:))*fs;
eeg(3,:) = ifft(eeg(3,:));
eeg(4,:) = ifftshift(EEG(4,:))*fs;
eeg(4,:) = ifft(eeg(4,:));
eeg(5,:) = ifftshift(EEG(5,:))*fs;
eeg(5,:) = ifft(eeg(5,:));

%plotting non-ideal filtered eeg signals
subplot(3,2,1)
plot(t, eeg(1,:))
xlabel('Time [s]')
ylabel('Amplitude')
title('Filtered EEG signal 1 - Time domain')
subplot(3,2,2)
plot(t, eeg(2,:))
xlabel('Time [s]')
ylabel('Amplitude')
title('Filtered EEG signal 2 - Time domain')
subplot(3,2,3)
plot(t, eeg(3,:))
xlabel('Time [s]')
ylabel('Amplitude')
title('Filtered EEG signal 3 - Time domain')
subplot(3,2,4)
plot(t, eeg(4,:))
xlabel('Time [s]')
ylabel('Amplitude')
title('Filtered EEG signal 4 - Time domain')
subplot(3,2,5)
plot(t, eeg(5,:))
xlabel('Time [s]')
ylabel('Amplitude')
title('Filtered EEG signal 5 - Time domain')

%% Question 2.9 equivilance with convoloution

%impulse response
impresp = h_f;

%convoloution of signalS
Conveq = conv(impresp, xdm(1,:))/fs;
%removing extra samples
Conveq = Conveq(1:length(xdm));

figure (6)
subplot(2,2,1)
plot (t, eeg(1,:))
xlabel('Time [s]')
ylabel('Amplitude')
title('EEG signal 1 - multiplied in frequency domain')
subplot(2,2,2)
plot(t, Conveq)
xlabel('Time [s]')
ylabel('Amplitude')
title('EEG signal 1 - Convoluted in time domain ')

%% Question 2.10 comparing signals

% comparison of the plots show the eeg signals are
% still noisey and need further filtering


%% Question 2.11 removing excess noise from signals

% by examining all signals, we can saftely determine band width between -60
% & 60 Hz

filtsig = zeros(5, length(EEG));


idealfilt = zeros(1, length(filtsig));
idealfilt(-60 < k & k < 60) = 1;
idealfilt(-56 < k & k < -54) = 0;
idealfilt(56 < k & k < 54) = 0;


%filtering EEG signals 
filtsig(1,:) = idealfilt.*EEG(1,:);
filtsig(2,:) = idealfilt.*EEG(2,:);
filtsig(3,:) = idealfilt.*EEG(3,:);
filtsig(4,:) = idealfilt.*EEG(4,:);
filtsig(5,:) = idealfilt.*EEG(5,:);


%bring EEG signals to time domain
filtsig(1,:) = ifftshift(filtsig(1,:))*fs;
filtsig(1,:) = ifft(filtsig(1,:));
filtsig(2,:) = ifftshift(filtsig(2,:))*fs;
filtsig(2,:) = ifft(filtsig(2,:));
filtsig(3,:) = ifftshift(filtsig(3,:))*fs;
filtsig(3,:) = ifft(filtsig(3,:));
filtsig(4,:) = ifftshift(filtsig(4,:))*fs;
filtsig(4,:) = ifft(filtsig(4,:));
filtsig(5,:) = ifftshift(filtsig(5,:))*fs;
filtsig(5,:) = ifft(filtsig(5,:));

%plotting against t vector
subplot(3,2,1)
plot(t, filtsig(1,:))
xlabel('Time [s]')
ylabel('Amplitude')
title('EEG signal - Astronaut 1')
subplot(3,2,2)
plot(t, filtsig(2,:))
xlabel('Time [s]')
ylabel('Amplitude')
title('EEG signal - Astronaut 2')
subplot(3,2,3)
plot(t, filtsig(3,:))
xlabel('Time [s]')
ylabel('Amplitude')
title('EEG signal - Astronaut 3')
subplot(3,2,4)
plot(t, filtsig(4,:))
xlabel('Time [s]')
ylabel('Amplitude')
title('EEG signal - Astronaut 4')
subplot(3,2,5)
plot(t, filtsig(5,:))
xlabel('Time [s]')
ylabel('Amplitude')
title('EEG signal - Astronaut 5')

%% Question 2.12 %% observing plots

%astronauts 1,3 and 4 are stable
%astronaut 5 is having a hard time
%astronaut 2 is on the verge of a mental breakdown








