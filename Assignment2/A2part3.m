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
first_image_data = data.sig(1,:);
imshow(reshape(first_image_data, 480, 640));
title('Corrupted image1');
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
ylabel('Magnitdue');
title('3 seconds of the received signal in Time Domain');

subplot(1,3,2)
image_shift = fftshift(fft(first_image_data))/fs;
plot(k, abs(image_shift));
SIG = fft(first_image_data);
xlabel('Frequency [Hz]');
ylabel('Magnitude');
title('Received signal in Frequency Domain');

subplot(1,3,3)
plot(k,abs(image_shift))
xlim([250 270])
xlabel('Frequency [Hz]');
ylabel('Magnitude');
title('Zoomed in Frequency Domain');

%% Question 3.4
 
T = candidateT(1,5);

Noisesig = estimateNoise(sig(1,:), T);
j = t(1:2*length(Noisesig));

figure(3)
%plot(j, Noisesig, 'r');
y = repmat(Noisesig.', [1, 2]);
plot(j, y,'r');
hold on

plot(t, first_image_data, 'b.');
xlim([0 3])
xlabel('Time [s]');
ylabel('Magnitude');
title('Image Signal and Periodic Noise Signal - 3 sec');
hold off
%% Question 3.5 
T = 1.469;
f = 1/T; % Fundamental frequency
a0 = (1/T).*sum(Noisesig'.*Ts); 
j = t(1:length(Noisesig));

% Using Fourier Series
FTSignal = a0;
N = 6; 
for n = 1:N
    an = (2/T).*sum(Noisesig'.*cos(2.*pi.*f.*n.*j))*Ts;
    bn = (2/T).*sum(Noisesig'.*sin(2.*pi.*f.*n.*j))*Ts;
    FTSignal = FTSignal + an.*cos(2.*pi.*f.*n.*j) + bn.*sin(2.*pi.*f.*n.*j);
end
figure(4)
plot(t(1:length(FTSignal)),FTSignal);
xlabel('Time [s]');
ylabel('Magnitude');
title('Fourier Transform Signal - 1.5 sec');

%% Question 3.6
% a0 represents the mean therefore that will be zero;
a0 = 0;
%% Question 3.7
%store in Noisesig_fs
Noisesig_fs = 0; %initialise Noisesig_fs
N = 20; % Increased number of harmonics

for n = 1:N
    an = (2/T).*sum(Noisesig'.*cos(2.*pi.*f.*n.*j))*Ts;
    bn = (2/T).*sum(Noisesig'.*sin(2.*pi.*f.*n.*j))*Ts;
    Noisesig_fs = Noisesig_fs + an.*cos(2.*pi.*f.*n.*j) + bn.*sin(2.*pi.*f.*n.*j);
end
% To calculate how many time the signal was repeated 
Rep = ceil(length(t) / 1469);
Noisesig_fs = repmat(Noisesig_fs, [1, Rep]);
Noisesig_fs = Noisesig_fs(1: length(t));
%% Question 3.8
figure(5)
plot(t(t<=j(end)), Noisesig_fs(t<=j(end)));
hold on
plot(j, Noisesig);
xlim([0 1.5])
xlabel('Time [s]');
ylabel('Magnitude');
title('Signal 1 and Modelled with 20 as the value of harmonics');
hold off
% Initial number of harmonics onsufficient to get accurate signal therefroe
% it is increased to 20
%% Question 3.9
%Site one with the periodic noise removed
figure(6)
De_noise = first_image_data - Noisesig_fs;
imshow(reshape(De_noise(1,:), 480, 640));
%plot(t, De_noise)
%xlabel('Time [s]');
%ylabel('Magnitude');
title('Site without periodic noise');
%comment on differences between this plot and first_image_data plot


%% Question 3.10

% Using Fourier transform to change the signal to frequency domain.
% Removing the bandwidth noise calculated and comverting back to time
% domain
De_noisefreq = fft(De_noise); 
De_noisefreq = fftshift(De_noisefreq)/fs;

sigfilt = zeros(1,length(De_noisefreq));

idealfilt = ones(1, length(sigfilt));
idealfilt(-270 < k & k <-250) = 0;
idealfilt(250 < k & k < 270) = 0;

sigfilt = idealfilt.*De_noisefreq;

figure(7)
im1 = ifftshift(sigfilt)*fs;
im1 = ifft(im1);

%Displaying the landing site 1 clearly
imshow(reshape(im1, 480, 640));
title('Clear image of site 1');

%% Question 3.11

% Performing the clearing for all four signals by creating custom filters
% for each signal
im2 = zeros(4, length(sig));

im2(1,:) = im1;
im2(2,:) = sig(2,:) - Noisesig_fs;
im2(3,:) = sig(3,:) - Noisesig_fs;
im2(4,:) = sig(4,:) - Noisesig_fs;

% Performing Fourier Transform
im1f = zeros(3,length(im2));
im1f(1,:) = fft(im2(2,:));
im1f(1,:) = fftshift(im1f(1,:))/fs; 
im1f(2,:) = fft(im2(3,:));
im1f(2,:) = fftshift(im1f(2,:))/fs; 
im1f(3,:) = fft(im2(4,:));
im1f(3,:) = fftshift(im1f(3,:))/fs; 

% Removing the custom bandwidth
finalfilt1 = ones(1, length(sigfilt));
finalfilt1(-245 < k & k < -225) = 0;
finalfilt1(225 < k & k < 245) = 0;

finalfilt2 = ones(1, length(sigfilt));
finalfilt2(-240 < k & k < -220) = 0;
finalfilt2(220 < k & k < 240) = 0;

finalfilt3 = ones(1, length(sigfilt));
finalfilt3(-255 < k & k < -235) = 0;
finalfilt3(235 < k & k < 255) = 0;


filt_im1 = zeros(3,length(im1f));

filt_im1(1,:) = im1f(1,:).*finalfilt1;
filt_im1(2,:) = im1f(2,:).*finalfilt2;
filt_im1(3,:) = im1f(3,:).*finalfilt3;


%Chaning back to time domain
im1t(1,:) = ifftshift(filt_im1(1,:))*fs;
im2(2,:) = ifft(im1t(1,:));

im1t(2,:) = ifftshift(filt_im1(2,:))*fs;
im2(3,:) = ifft(im1t(2,:));

im1t(3,:) = ifftshift(filt_im1(3,:))*fs;
im2(4,:) = ifft(im1t(3,:));

%Displaying the landing sites
figure(8)
subplot(2,2,1)
imshow(reshape(im2(1,:), 480, 640));
title('Site 1');
subplot(2,2,2)
imshow(reshape(im2(2,:), 480, 640));
title('Site 2');
subplot(2,2,3)
imshow(reshape(im2(3,:), 480, 640));
title('Site 3');
subplot(2,2,4)
imshow(reshape(im2(4,:), 480, 640));
title('Site 4');

%% Question 3.12

% Displaying the images seperately so the coordinates can be seen clearly
figure(9)
imshow(reshape(im2(1,:), 480, 640));
title('Site 1');

figure(10)
imshow(reshape(im2(2,:), 480, 640));
title('Site 2');

figure(11)
imshow(reshape(im2(3,:), 480, 640));
title('Site 3');

figure(12)
imshow(reshape(im2(4,:), 480, 640));
title('Site 4');
    




