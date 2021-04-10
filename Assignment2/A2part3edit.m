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
%% Question 3.1
data = matfile('A2P3Data.mat');
figure(1)
title('Corrupted image1');
first_image_data = data.sig(1,:);
imshow(reshape(first_image_data, 480, 640));

%% Question 3.2
samples = length(sig);
fs = 1000;
Ts = 1/1000;
% Creating the time vector
t = linspace(0, Ts * length(first_image_data), length(first_image_data) + 1); t(end) = []; 
% Creating the Frequency vector
k = linspace(-fs/2,fs/2, length(first_image_data) + 1); k(end) = []; 

%% Question 3.3
figure(2)
subplot(1,3,1)
plot(t, first_image_data);
xlim([0 3]);
xlabel('Time [s]');
ylabel('Amplitude');
title('Time Domain');

subplot(1,3,2)
image_shift = fftshift(fft(first_image_data))/fs;
plot(k, abs(image_shift));
SIG = fft(first_image_data);

subplot(1,3,3)
plot(k,abs(image_shift))
xlim([250 270])

% The periodic signal can be observed in the time frame of 0.292 to 1.76
% seconds giving us a period of 1.469 seconds


%% Question 3.4
% EXPLAIN WHY U PICKED THIS
T = candidateT(1,5);


Noisesig = estimateNoise(sig(1,:), T);
%Noisesig = repmat(Noisesig, [3 1]);
j = t(1:length(Noisesig));

% Repeat this viz for each period in CandidateT
figure(3)
plot(j, Noisesig, 'r');
hold on
plot(t, first_image_data, 'b.');
xlim([0 3])
hold off

%% Question 3.5 Trig
T = 1.469; % From part 3.4
f = 1/T; % Fundamental frequency
a0 = (1/T).*sum(Noisesig'.*Ts); 
j = t(1:length(Noisesig));
%tStep = period(2) - period(1);

FTSignal = a0;
N = 6; 
for n = 1:N
    an = (2/T).*sum(Noisesig'.*cos(2.*pi.*f.*n.*j))*Ts;
    bn = (2/T).*sum(Noisesig'.*sin(2.*pi.*f.*n.*j))*Ts;
    FTSignal = FTSignal + an.*cos(2.*pi.*f.*n.*j) + bn.*sin(2.*pi.*f.*n.*j);
end

figure(4)
plot(t(1:length(FTSignal)),FTSignal);
%% Question 3.6
% Can only be done once the other questions are complete 
% a0 represents the mean therefore that will be zero;
a0 =0;

%% Question 3.7
%store in Noisesig_fs
Noisesig_fs = 0; %initialise Noisesig_fs
N = 20;
for n = 1:N
    an = (2/T).*sum(Noisesig'.*cos(2.*pi.*f.*n.*j))*Ts;
    bn = (2/T).*sum(Noisesig'.*sin(2.*pi.*f.*n.*j))*Ts;
    Noisesig_fs = Noisesig_fs + an.*cos(2.*pi.*f.*n.*j) + bn.*sin(2.*pi.*f.*n.*j);
end

Rep = ceil(length(t) / 1469);
Noisesig_fs = repmat(Noisesig_fs, [1, Rep]);
Noisesig_fs = Noisesig_fs(1: length(t));

%% Question 3.8

figure(5)
plot(t(t<=j(end)), Noisesig_fs(t<=j(end)));
hold on
plot(j, Noisesig);
xlim([0 1.5])
hold off

% Initial number of harmonics insufficient to get accurate signal therefore
% it is increased to 20

%% Question 3.9
figure(6)

De_noise = first_image_data - Noisesig_fs;
imshow(reshape(De_noise(1,:), 480, 640));

%% Question 3.10
% Don't really know how to do this

%% Question 3.11
% Can only be done once the other questions are complete 

%% Question 3.12
% Can only be done once the other questions are complete 





