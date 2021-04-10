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
xlabel('Frequency [Hz]');
ylabel('Magnitude');
title('Frequency Domain');

SIG = fft(first_image_data);

subplot(1,3,3)
plot(k,abs(image_shift))
xlim([250 270])
xlabel('Frequency [Hz]');
ylabel('Magnitude');
title('Bandlimited Noise');

% The periodic signal can be observed in the time frame of 0.292 to 1.76
% seconds giving us a period of 1.469 seconds


%% Question 3.4
T = candidateT(1,5);

Noisesig = estimateNoise(sig(1,:), T);
Noisesig = repmat(Noisesig, [3 1]);
j = t(1:length(Noisesig));

figure(3)
plot(j, Noisesig, 'r');
hold on
plot(t, first_image_data, 'b.');
xlim([0 3])
hold off
xlabel('Time [s]');
ylabel('Amplitude');
title('Image Signal and Periodic Noise Signal - 3 sec');
legend('Periodic Noise Signal','1st Received Image Signal');

%% Question 3.5
% Decided to use Trigonometric Fourier Series 
T = 1.469; % From part 3.4
f = 1/T; % Fundamental frequency
a0 = (1/T).*sum(Noisesig'.*Ts); 
j = t(1:length(Noisesig));
%tStep = period(2) - period(1);

FTSignal = a0;
N = 100; 
for n = 1:N
    an = (2/T).*sum(Noisesig'.*cos(2.*pi.*f.*n.*j))*Ts;
    bn = (2/T).*sum(Noisesig'.*sin(2.*pi.*f.*n.*j))*Ts;
    FTSignal = FTSignal + an.*cos(2.*pi.*f.*n.*j) + bn.*sin(2.*pi.*f.*n.*j);
end
figure
plot(t(1:length(FTSignal)),FTSignal);
hold on
plot(j, Noisesig, 'r');
hold off
xlabel('Time [s]');
ylabel('Amplitude');
title('Image Signal and Periodic Noise Signal - 3 sec');
legend('Fourier Approximation - 6 Harmonics','Periodic Noise Signal');