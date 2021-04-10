%% Assignment 2, Part 3 (Choosing a landing site)
%  Do not change before line 32
%  You will need to have generated A2P3Data.mat from 
%  GenerateAssignment2Data.m before working with this file.

%  Clearing and preparing the workspace
clear; clc; close all;

%  Load assignment data from A2P3Data.mat.
load('A2P3Data.mat');  

%=================================================================%
%
% Refer to the assignment sheet for details.
% Names of the variables are important,
% e.g. 'a1' is considered a different variable to 'A1'.
% Make sure variables have been declared as they appear in the brief.
%
% t - time domain vector
% k - frequency vector
% SIG - Fourier Transform
% T - selected value
% Noisesig - estimated noise signal
% a0,an,bn - Trigonometric Fourier Series Coefficients
% OR
% c0,cn - Complex Fourier Series Coefficients
% Noisesig_fs - approximation of noise
% im1 - image 1
% im2 - image 2
% 
%====Enter your code below this line================================
%% 3.1 Viewing the Noisy Image
figure(1)
imshow(reshape(sig(1,:), 480, 640)) %show noisy image
title('The Noisy Image')

%% 3.2 Creating vectors for time and frequency 
samples = length(sig);
fs = 1000; %provided value for sampling frequency
Ts = 1/fs; %sampling period (period is inversly proportional to frequency
secs = samples/fs; %total time in seconds of the signal 
t = linspace(0, secs, samples+1); t(end) = []; %time vector
k = linspace(-fs/2,fs/2, samples+1); k(end) = []; %Frequency vector

%% 3.3

% Plots the Received Signal in the time and freqency domains
figure(2)
plot(t, sig(1,:));
title('The First Received Signal over Three Seconds')
xlabel('Time (s)')
ylabel('Amplitude')
xlim([0 3])

SIG = fft(sig(1,:));
SIG_shift = fftshift(abs(SIG)/fs);
figure(3)
plot(k, SIG_shift);
title('The First Received Signal in the Frequency Domain')
xlabel('Frequency (Hz)')
ylabel('Magnitude')

%% 3.4 

%Calculates an estimate of the interfering periodic noise given the vector candidateT and plots against the received signal
for i = 1:length(candidateT) 
    Noisesig1 = estimateNoise(sig(1,:),candidateT(i));
    Noisesig = repmat(Noisesig1, [3 1]); %3 seconds

    figure(4)
    subplot(5,2,i)
    plot(t,sig(1,:),'b') 
    hold on 
    a = t(1:length(Noisesig));
    plot(a,Noisesig,'r') %plots the noise profile against corresponding time vector a
    xlim([0,3])
    title(sprintf('Periodic Noise Signal %d Vs time',i)) 
    % (and the Recieved signal)
    xlabel('Time (s)')
    ylabel('Magnitude')
    
end
h = legend('Signal','Estimated Noise');
set(h, 'Position', [0.06 0.9 0 0])
 


%% 3.5

% Takes the best candidate period value and uses that value to calculate the Trigonometric Fourier coefficients

% Given in previous section, the noisesignal fits best is with
% candidateT value 8, which is 1109 
T = candidateT(8);

nSig = estimateNoise(sig(1,:),T);
tn = t(1:length(nSig));
T = T/fs;
nSig = nSig';

% TRIG COEFFICIENTS
Nt = 6;
nt = 1:Nt; %1:1:Nt
n = (1:Nt).';
f0 = 1/T;

a0 = (1/T) * sum(nSig)/fs;
an = (2/T) * nSig * cos(2*pi*f0*n*tn).'/fs;
bn = (2/T) * nSig * sin(2*pi*f0*n*tn).'/fs;

%% 3.6
% Removes the DC offset by setting a0 to 0, therefore creating the mean 
    % of the periodic noise signal equal to 0
a0 = a0*0; % set to zero

%% 3.7 

% Calculates the Trigonometric Fourier Series with 6 harmonics
w0 = (2*pi)/T;
Noisesig_fs = 0; %initialising Noisesig_fs before the for loop
for i = 1:Nt
    harm = a0 + an(i)*cos(i*w0*t)+bn(i)*sin(i*w0*t); %the trig fourier series excluding a0
    Noisesig_fs = Noisesig_fs + harm; 
end

%% 3.8 

% Plots the Fourier Series approximation and the received signal in the time domain
figure(5)
subplot(2,1,1)
plot(t,Noisesig_fs)
title('TFS Approximation with 6 Harmonics');
xlabel('Time (s)');
ylabel('Amplitude'); 
xlim([0 T])
ylim([-6 12])

subplot(2,1,2)
plot(t, sig(1,:));
title('The First Received Signal over the First Period')
xlabel('Time (s)')
ylabel('Amplitude')
xlim([0 T])
ylim([-6 12])

%% 3.9

% Creates the noisesig variable and replicates it to be a 1x307200 vector
repetitions = round(secs/T); %total time of signal over period 
Noisesig_fsR = repmat(Noisesig_fs(1,:),1,repetitions);
Noisesig_fs1 = Noisesig_fsR(1:length(t)); 

% removes periodic noise signal from received signal
im1 = sig(1,:) - Noisesig_fs1;
figure(6) 
imshow(reshape(im1(1,:),480,640)); %shows denoised image 1
title("Recovered Image after Approximate Noise Removed")

%frequency spectrum
IM1 = fft(im1(1,:));
IM1_shift = fftshift(abs(IM1))/fs;

% Plots the denoised and received signals in the frequency domain
figure(7)
subplot(2,1,1)
plot(k,abs(fftshift(SIG)/fs))
title('Received Signal Frequency Spectrum')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
subplot(2,1,2)
plot(k,IM1_shift)
title('Frequency Spectrum of Denoised Image')
xlabel('Frequency (Hz)')
ylabel('Amplitude')

% This whole section calculates trig coefficients for harmonics of
    % 10, 20, 30, 50, and 100, as well as calculating the denoised signal

%10 Coefficients
[a0_10, an_10, bn_10, Noisesig_fs_10] = TrigFSM(nSig, 10, T, fs, tn, t);

Noisesig_fsR_10 = repmat(Noisesig_fs_10(1,:),1,repetitions);
Noisesig_fs1_10 = Noisesig_fsR_10(1:length(t)); 

im1_10 = sig(1,:) - Noisesig_fs1_10;
%frequency spectrum
IM1_10 = fft(im1_10(1,:));
IM1_10_shift = fftshift(abs(IM1_10))/fs;

%20 Coefficients
[a0_20, an_20, bn_20, Noisesig_fs_20] = TrigFSM(nSig, 20, T, fs, tn, t);

Noisesig_fsR_20 = repmat(Noisesig_fs_20(1,:),1,repetitions);
Noisesig_fs1_20 = Noisesig_fsR_20(1:length(t));

im1_20 = sig(1,:) - Noisesig_fs1_20;
%frequency spectrum
IM1_20 = fft(im1_20(1,:));
IM1_20_shift = fftshift(abs(IM1_20))/fs;

%30 Coefficients
[a0_30, an_30, bn_30, Noisesig_fs_30] = TrigFSM(nSig, 30, T, fs, tn, t);

Noisesig_fsR_30 = repmat(Noisesig_fs_30(1,:),1,repetitions);
Noisesig_fs1_30 = Noisesig_fsR_30(1:length(t)); 

im1_30 = sig(1,:) - Noisesig_fs1_30;
%frequency spectrum
IM1_30 = fft(im1_30(1,:));
IM1_30_shift = fftshift(abs(IM1_30))/fs;

%50 Coefficients
[a0_50, an_50, bn_50, Noisesig_fs_50] = TrigFSM(nSig, 50, T, fs, tn, t);

Noisesig_fsR_50 = repmat(Noisesig_fs_50(1,:),1,repetitions);
Noisesig_fs1_50 = Noisesig_fsR_50(1:length(t)); 

im1_50 = sig(1,:) - Noisesig_fs1_50;
%frequency spectrum
IM1_50 = fft(im1_50(1,:));
IM1_50_shift = fftshift(abs(IM1_50))/fs;

%100 Coefficients
[a0_100, an_100, bn_100, Noisesig_fs_100] = TrigFSM(nSig, 100, T, fs, tn, t);

Noisesig_fsR_100 = repmat(Noisesig_fs_100(1,:),1,repetitions);
Noisesig_fs1_100 = Noisesig_fsR_100(1:length(t)); 

im1_100 = sig(1,:) - Noisesig_fs1_100;
%frequency spectrum
IM1_100 = fft(im1_100(1,:));
IM1_100_shift = fftshift(abs(IM1_100))/fs;


% Plots the frequency spectra of all denoised signals calculated above
figure(8)
sgtitle("Bandlimited Noise Removed");
subplot(2,3,1)
plot(k, abs(fftshift(IM1)/fs))
title('IM1 6 Harmonics')
set(gca,'XTick',(-500:100:500))

subplot(2,3,2)
plot(k, abs(fftshift(IM1_10)/fs))
title('IM1 10 Harmonics')
set(gca,'XTick',(-500:100:500))

subplot(2,3,3)
plot(k, abs(fftshift(IM1_20)/fs))
title('IM1 20 Harmonics')
set(gca,'XTick',(-500:100:500))

subplot(2,3,4)
plot(k, abs(fftshift(IM1_30)/fs))
title('IM1 30 Harmonics')
set(gca,'XTick',(-500:100:500))

subplot(2,3,5)
plot(k, abs(fftshift(IM1_50)/fs))
title('IM1 50 Harmonics')
set(gca,'XTick',(-500:100:500))

subplot(2,3,6)
plot(k, abs(fftshift(IM1_100)/fs))
title('IM1 100 Harmonics')
set(gca,'XTick',(-500:100:500))


% Plots the denoised images
figure(9) %Denoised Image
sgtitle(["Denoised Image im1",""]);
subplot(2,3,1)
imshow(reshape(im1(1,:),480,640));
title("im1 6 Harmonics") 

subplot(2,3,2)
imshow(reshape(im1_10(1,:),480,640));
title("im1 10 Harmonics")

subplot(2,3,3)
imshow(reshape(im1_20(1,:),480,640));
title("im1 20 Harmonics")

subplot(2,3,4)
imshow(reshape(im1_30(1,:),480,640));
title("im1 30 Harmonics")

subplot(2,3,5)
imshow(reshape(im1_50(1,:),480,640));
title("im1 50 Harmonics")

subplot(2,3,6)
imshow(reshape(im1_100(1,:),480,640));
title("Denoised Image with 100 Harmonics")

% Plots the Trig Fourier series approximations for 6 and 50 harmonics 
    %and the received signal in the time domain
figure(22)
sgtitle("Trigonometric Fourier Approximation");
subplot(3,1,1)
plot(t,Noisesig_fs)
title('TFS Approximation with 6 Harmonics');
xlabel('Time (s)');
ylabel('Amplitude'); 
xlim([0 T])
ylim([-6 12])

subplot(3,1,2)
plot(t, Noisesig_fs_50)
title('TFS Approximation with 50 Harmonics');
xlabel('Time (s)');
ylabel('Amplitude'); 
xlim([0 T])
ylim([-6 12])

subplot(3,1,3)
plot(t, sig(1,:));
title('The First Received Signal over the First Period')
xlabel('Time (s)')
ylabel('Amplitude')
xlim([0 T])
ylim([-6 12])

%% 3.10 

% Calculates the denoised signal for all harmonic values used
im1_6(1,:) = sig(1,:) - Noisesig_fs1;
im1_10(1,:) = sig(1,:) - Noisesig_fs1_10;
im1_20(1,:) = sig(1,:) - Noisesig_fs1_20;
im1_30(1,:) = sig(1,:) - Noisesig_fs1_30;
im1_50(1,:) = sig(1,:) - Noisesig_fs1_50;
im1_100(1,:) = sig(1,:) - Noisesig_fs1_100;


% Applies a bandstop filter betweeh 100 and 270Hz in the positive and negative frequencies
im2_6 = bandstop(im1_6,[100 270],fs);
im2_10 = bandstop(im1_10,[100 270],fs);
im2_20 = bandstop(im1_20,[100 270],fs);
im2_30 = bandstop(im1_30,[100 270],fs);
im2_50 = bandstop(im1_50,[100 270],fs);
im2_100 = bandstop(im1_100,[100 270],fs);

% This section plots the denoised images after the filter is applied
figure(10) 
imshow(reshape(abs(im2_6(1,:)),480,640))
title('6 Harmonics')

figure(11)
imshow(reshape(abs(im2_10(1,:)),480,640))
title('10 Harmonics')

figure(12)
imshow(reshape(abs(im2_20(1,:)),480,640))
title('20 Harmonics')

figure(13)
imshow(reshape(abs(im2_30(1,:)),480,640))
title('30 Harmonics')

figure(14)
imshow(reshape(abs(im2_50(1,:)),480,640))
title('50 Harmonics')

figure(15)
imshow(reshape(abs(im2_100(1,:)),480,640))
title('100 Harmonics')
%% 3.11

%This applies everything learned in the above sections for all 4 received signals
for i = 1:4
   im1(i,:) = sig(i,:) - Noisesig_fs1_50;
   IM1(i,:) = fft(im1(1,:));
   IM1_shift(i,:) = fftshift(abs(IM1(i,:)))/fs;
   
   figure(16)
   subplot(4,1,i)
   plot(k,IM1_shift(i,:))
   title(sprintf('Frequency Spectrum of Denoised Image %d',i))
   xlabel('Frequency (Hz)')
   ylabel('Amplitude')
   
   % Plots each image in one plane
   im2 = zeros(4,307200);
   im2(i,:) = bandstop(im1(i,:),[100 270],fs);
   figure(17)
   subplot(2,2,i)
   imshow(reshape(abs(im2(i,:)),480,640))
   title(sprintf('Landing Site %d',i))
   
   %3.12
   % Plots each image separately to get a better look at the images and see the navigation markers
   figure(17+i)
   imshow(reshape(abs(im2(i,:)),480,640))
   title(sprintf('Landing Site %d',i))  
end
%% Functions

% Function to calculate the Trig Coefficients and Fourier Series for part 3.9
function [a0f, anf, bnf, Noisesig_fsf] = TrigFSM(sig, N, T, fs, t1, t)
    % TRIG COEFFICIENTS
    nf = (1:N).';
    f0f = 1/T;

    a0f = (1/T) * sum(sig)/fs;
    anf = (2/T) * sig * cos(2*pi*f0f*nf*t1).'/fs;
    bnf = (2/T) * sig * sin(2*pi*f0f*nf*t1).'/fs;

    w0 = (2*pi)/T;
    Noisesig_fsf = 0; %initialising Noisesig_fs before the for loop
    for i = 1:N
        harmf = anf(i)*cos(i*w0*t)+bnf(i)*sin(i*w0*t);
        Noisesig_fsf = Noisesig_fsf + harmf; 
    end
end