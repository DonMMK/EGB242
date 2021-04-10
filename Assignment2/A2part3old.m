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
first_image_data = data.sig(1,:);
imshow(reshape(first_image_data, 480, 640));

%% Question 3.2
fs = 1000;
Ts = 1/1000;
t = linspace(0, Ts * length(first_image_data), length(first_image_data) + 1); t(end) = []; % time vector  
k = linspace(-fs/2,fs/2, length(first_image_data) + 1); k(end) = []; % frequency vector

%% Question 3.3
figure
subplot(2,1,1)
plot(t, first_image_data);
xlim([0 3]);
subplot(2,1,2)
image_shift = fftshift(fft(first_image_data))/fs;
plot(k, abs(image_shift));
SIG = fft(first_image_data);
% The plots look different from the ones showed in the assignment walkthrough 
%I must of done something wrong here
% Second part of the question is for the report I have a general idea but
% not too sure

%% Question 3.4
T = (1630)*10^-3; % 1630 is the approximate period that 
%I got based on the plots that I have generated 
Noisesig = estimateNoise(first_image_data, T*fs);
% Did this based on the assignment walkthrough but comes out as an error 
figure
plot(first_image_data, Noisesig);

%% Question 3.5
% Decided to use Complex Fourier Series Coefficients
N = 6;
period = linspace(0,6, Noisesig + 1); period(end) = [];
tStep = period(2) - period(1);
T = 6;

n = (-N:N)'; n(n == 0) = [];

    subfunction = 2 * pi * n * period * (1/T);
    
    % Finding Coefficients
    c0 = 1/T * sum(Noisesig) * tStep;
    cn = 1/T * Noisesig * exp(1j*-subfunction).'*tStep;
% Haven't listed the coefficients and explained the process as I can only
% do that once the code is correct
%% Question 3.6
% Can only be done once the other questions are complete 

%% Question 3.7
FS = cn * exp(1j*subfunc);
tPeriod = linspace(0, Ts, 307200 + 1); tPeriod(end) = [];
Noisesig_fs = [(FS) (FS)];
figure
plot(tPeriod, Noisesig_fs);
%Not sure if I am right with this
%% Question 3.8
% Can only be done once the other questions are complete and doesn't
% require any code

%% Question 3.9
De_noise = real(noiseSound - Noisesig_fs');
im1 = [De_noise];

%% Question 3.10
% Don't really know how to do this

%% Question 3.11
% Can only be done once the other questions are complete 

%% Question 3.12
% Can only be done once the other questions are complete 





