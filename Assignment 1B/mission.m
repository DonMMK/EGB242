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
%   Ts              Sampling period
%   t               time vector
%   MUX             Fourier transform of muxSignal
%   MUX_shift       Shifted and scaled Fourier Transform
%   f               Frequency vector
%   freqshift       Frequency shifts (vector)
%   MagSpec         amplitude of peaks (vector)
%   PhaseSpec       phases of peaks (vector)
%   xdm             output of FDMDemux (matrix)
%   XDM             Fourier transform of xdm (matrix)
%   B               Bandwidth
%   filteredoutput  Filtered Output Signal
%   decodedtext     The Decoded Text Output
% ---------------------------------------------
%====Enter your code below this line================================

%% Section B2.1

% Determine the sampling Period
Ts = 1/fs;

% Construct the time domain vector 
t = linspace(0, Ts*length(muxSignal), length(muxSignal)+1 );t(end) = [];

% Plotting muxSignal against Time
figure(1);
plot(t, muxSignal);
grid on
xlabel('Time (s)');
ylabel('Magnitude');
title('Plotting muxSignal');

%% Section B2.2

% Computing the Fourier Transform
MUX = fft(muxSignal);

% Constructing an appropriate frequency vector
f = linspace(-fs/2 , fs/2 , length(muxSignal) + 1); f(end) = [];

% Generating MUX shift
MUX_shift = fftshift(MUX);

% Plotting the graph of frequency against magnitude
figure(2);
plot(f, abs(MUX_shift)/fs);
grid on
title('Plot of MUX magnitude spectrum')
xlabel('Frequency (Hz)')
ylabel('Magnitude')

%% Section B2.3

%store in freqshift
% Locating and storing the locations of the peaks
[Peaks, PeakLocations] = findpeaks(abs(MUX_shift)/fs, 'MinPeakHeight', 0.4);

% To get only positive values and sorting in the correct order
index = sort(PeakLocations(6:10)); 

% Getting the required frequency shifted values
freqshift = f(index);

%% Section B2.4

% Find magnitude and store in MagSpec and find phase and save in PhaseSpec
% The values in these vectors need to correspond to their respective frequencies in freqshift.
MagSpec = abs(MUX_shift(index))/fs;
PhaseSpec = angle(MUX_shift(index)); 

% Plot locations and magnitudes of the frequency shifts using red circles on the magnitude spectrum of MUX
figure(3);
hold on
plot(f, abs(MUX_shift)/fs);
plot(freqshift , MagSpec , 'ro')
grid on
xlabel('Frequency (Hz)')
ylabel('Magnitude');
title('Plot of magnitude against frequency');
hold off
%% Section B2.5

% FDMDemux function implemented into MATLAB so that mission.m can run
% properly
    xdm = zeros(length(freqshift),length(t));
    
    for i = 1:length(freqshift)
        xdm(i,:) = muxSignal.*cos(2*pi*freqshift(i)*t+PhaseSpec(i))*MagSpec(i); 
    end

% frequency shifting module called FDMDemux.m
[xdm] = FDMDemux(muxSignal,t,MagSpec,freqshift,PhaseSpec);

%% Section B2.6

% Compute the Fourier transform for each data stream in xdm (row by row)
% Store the result in the matrix XDM.
XDM = zeros( size(xdm));

 for count = 1: size(xdm , 1)
     XDM(count, :) = fft(xdm(count, : ));
     
 end   
     
% Plot the magnitude spectrum for each data stream. 
figure(4);
for count = 1: size(XDM)
   subplot(2, 3, count);
   plot(f, abs(fftshift(XDM(count, :)) / fs));
   title('Magnitudes from Data stream')
end
%% Section B2.7
 
% Storing the bandwidth value in the vector B (Bandwidth)
B = 2000;
% All possible values of B = [5200 2500 2200 2400 5400];

% Creating the required filter
Filter = zeros(1,length(t));
Filter( -B<f &  f < B) = 1;

Filtered_Signal1 = fftshift(Filter).*(XDM(1,:));
Filtered_Signal2 = fftshift(Filter).*(XDM(2,:));
Filtered_Signal3 = fftshift(Filter).*(XDM(3,:));
Filtered_Signal4 = fftshift(Filter).*(XDM(4,:));
Filtered_Signal5 = fftshift(Filter).*(XDM(5,:));


figure(5)
subplot(5,1,1)
plot(f, fftshift(abs(Filtered_Signal1))/fs);
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Filtered Data from stream 1');
subplot(5,1,2)
plot(f, fftshift(abs(Filtered_Signal2))/fs);
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Filtered Data from stream 2');
subplot(5,1,3)
plot(f, fftshift(abs(Filtered_Signal3))/fs);
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Filtered Data from stream 3');
subplot(5,1,4)
plot(f, fftshift(abs(Filtered_Signal4))/fs);
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Filtered Data from stream 4');
subplot(5,1,5)
plot(f, fftshift(abs(Filtered_Signal5))/fs);
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Filtered Data from stream 5');

figure(6)
subplot(5,1,1)
signal1 = ifft(Filtered_Signal1) - mean(ifft(Filtered_Signal1));
plot(t, signal1)
xlabel('Time (s)')
ylabel('Magnitude')
title(' Decoded signal stream 1')

subplot(5,1,2)
signal2 = ifft(Filtered_Signal2) - mean(ifft(Filtered_Signal2));
plot(t, signal2)
xlabel('Time (s)')
ylabel('Magnitude')
title(' Decoded signal stream 2')

subplot(5,1,3)
signal3 = ifft(Filtered_Signal3) - mean(ifft(Filtered_Signal3));
plot(t, signal3)
xlabel('Time (s)')
ylabel('Magnitude')
title(' Decoded signal stream 3')

subplot(5,1,4)
signal4 = ifft(Filtered_Signal4) - mean(ifft(Filtered_Signal4));
plot(t, signal4)
xlabel('Time (s)')
ylabel('Magnitude')
title(' Decoded signal stream 4')

subplot(5,1,5)
signal5 = ifft(Filtered_Signal5) - mean(ifft(Filtered_Signal5));
plot(t, signal5)
xlabel('Time (s)')
ylabel('Magnitude')
title(' Decoded signal stream 5')

%% Section B2.8

% Signal 1 
sound(signal1,fs)
% Signal 2
text2 = A1BTextdecode(signal2,fs)
% Signal 3
text3 = A1BTextdecode(signal3,fs)
% Signal 4
text4 = A1BTextdecode(signal4,fs)
% Signal 5
sound(signal5,fs)

decodedtext = [text2 text3 text4 ];
