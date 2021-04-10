%% Assignment 1 - Part B, Section B3 (Combining Signals)
%  Do not change before line 40
%  You will need to have generated Data1B.mat from 
%  GenerateDataAssignment1B.m before working with this file.

%  Clearing and preparing the workspace
clear; clc; close all;

%  Load assignment data from Data1B.mat
load('Data1B.mat', 'fs', 'muxSignal');% 

%  The variable loaded are:
%     fs      Sampling frequency for Section B3
%  muxSignal  Multiplexed music signals, for Section B3.

%==================================================================
%
% Names of variables you will need for marking.
% Refer to the assignment sheet for details.
% Names of the variables are important,
% e.g. 'a1' is considered a different variable to 'A1'.
% Make sure variables have been declared as they appear in the list below.
% ---------------------------------------------
% ===== Part 2 =====
%     Ts         Sampling period
%     t          time vector
%     MUXF        Fourier transform of muxSignal
%   MUXF_shift    Shifted and scaled Fourier Transform
%     k          Frequency vector
%   freqshifting Frequency shifts (vector)
%   MagSpec      mahnitude of peaks (vector)
%   PhaseSpec    phases of peaks (vector)
%     xdm        output of FDMDemux (matrix)
%     XDM        Fourier transform of xdm (matrix)
%     B          Bandwidth
% filteredoutput Filtered Output Signal
% decodedtext    The Decoded Text Output
% ---------------------------------------------
%====Enter your code below this line================================


%% Plotting The Signal Of The Raw Data Stream %% DONE %%%%%

%The Raw Data Stream Of Music And Text Signals Are Stored In 'muxSignal'
%The Sampling Frequency of 'muxSignal' Is Stored in 'fs' (96000Hz/96kHz)

%Sampling Period Of The Raw Data Stream
Ts = 1/fs;

%Creating An Appropriate Time Vector Using The Sampling Period
%t = linspace(-Ts/2, Ts/2, length(muxSignal) + 1); t(end) = [];
%t = linspace(0, 10, 1056000+1); t(end) = [];
t = linspace(0, Ts*length(muxSignal),length(muxSignal) + 1); t(end) =[];

%Plotting The Signal Containing The Raw Data Stream
figure(1);
plot(t, muxSignal);
grid on

%Labelling Axes And Adding Title 
xlabel('Time (s)');
ylabel('Magnitude');
title('Raw Data Stream (muxSignal) vs Time');


%% Fourier Transform And Magnitude Spectrum %% DONE %%%%%

%Computing The Fourier Transform Of The Raw Data Stream
MUXF = fft(muxSignal);

%Creating An Appropriate Frequency Vector For MUXF
k = linspace(-fs/2, fs/2, length(muxSignal) + 1); k(end) = [];

%Plotting The Fourier Transform Of The Raw Data Stream
figure(2);
plot(k, real(MUXF));
grid on

%Labelling Axes And Adding Title
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Fourier Transform Of Raw Data Stream (MUXF) vs Frequency');

%Computing The Shifted And Scaled Fourier Trasform Of muxSignal
MUXF_shift = abs(fftshift(MUXF/fs));

%Plotting The Magnitude Spectrum Of MUXF Against Frequency
figure(3);
plot(k, MUXF_shift);
grid on

%Labelling Axes And Adding Title
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Magnitude Spectrum Of MUXF vs Frequency');


%% Frequency Shifts   %% DONE %%%%%

%Identifying All Peaks(Magnitude) And Their Respective Locations(Frequency)
%Using Inbuilt Matlab 'findpeaks' Function
[peaks,locations] = findpeaks(MUXF_shift,k,"MinPeakHeight",0.4);

%Getting Posititve Frequencies Only
%Sorting The Values In Ascending Order
%Rounding Off To The Nearest Integer
%Storing The Frequencies In A Row Vector Named 'freqshifting'
freqshifting = round(sort(locations(locations>0)));


%% Magnitude Spectrum And Phase Spectrum %% DONE %%%%%

%Finding Magnitude Of The Peaks And Storing In MagSpec
MagSpec = peaks(locations>0);

%Finding Phase Of The Peaks And Storing In PhaseSpec
A = fftshift(MUXF);
PhaseSpec = angle(A(locations>0));

%Plotting Magnitude Spectrum of MUXF Against Frequency
figure(4); 
plot(k, MUXF_shift);
hold on

%Plotting An Overlay Of Location And Magnitude Of Frequency Shifts
plot(freqshifting,MagSpec,'ro'); %Red Circles
grid on

%Labelling Axes And Adding Title
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Magnitude Spectrum Of MUXF vs Frequency');


%% Frequency Shifting Module 'FDMDemux.m'

[xdm] = FDMDemux(muxSignal,t,MagSpec,freqshifting,PhaseSpec);


%% B2.6

XDM = zeros(size(xdm)); 
for i = 1:size(xdm,1)
    XDM(i,:) = fft(xdm(i,:)); 
end 

figure(5); 
for i = 1:size(XDM)
    subplot(2,3,i);
    plot(k, abs(fftshift(XDM(i,:))/fs)); 
    title(['Magnitude Spectrum At Data Set ',num2str(i)]); 
end


%% B2.7

B = [4912,2046,2386,5283];

filteredoutput = zeros(size(xdm));

for i =1:size(filteredoutput,1)
    filteredoutput(i,:) = A1BLPF(xdm(i,:),fs,(B(i)/2));
end

filteredoutput(1,:) = filteredoutput(1,:) + 0.02;
filteredoutput(2,:) = filteredoutput(2,:) + 0.01244;
filteredoutput(3,:) = filteredoutput(3,:) + 0.1147;
filteredoutput(4,:) = filteredoutput(4,:) - 0.04917;


%% B2.8

% Signal 1
text1 = A1BTextdecode(filteredoutput(2,:),Ts);

% % Signal 2
% sound(filteredoutput(1,:),fs)

% Signal 3
text2 = A1BTextdecode(filteredoutput(3,:),Ts);

% % Signal 4
% sound(filteredoutput(4,:),fs)

decodedtext = [text1 text2] % Decoded Text Output