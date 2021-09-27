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


%% --------------3.1 View the noisy image-----------------

%converting the signal into 2D and displaying the first noisy image 
sig1 = sig(1,:);
figure (1)
title('Signal 1 currupted image')
imshow(reshape(sig1,480,640));

%% --------------3.2 Reference vector -------------------

%Sampling Variables 
samples = length(sig);
fs = 1000;
Ts = 1/fs;

%Time vector

t=linspace(0,Ts*samples,samples+1);
t(end) =[];

%Frequency vector 

k = linspace(-fs/2,fs/2,samples+1);
k(end) =[];



%% --------------3.3 Visualise the received signal-------------


%First received signal in time domain (Fisrt 3 seconds)

figure(2)
plot(t,sig1,'k');
xlim([0 3]) 
xlabel('Time(Seconds)');
ylabel('Magnitude');
title(' 3 secods of the first received signal in time domain');

% performing Fourier transform on signal
SIG = fft(sig1);

%First received signal in frequency domain 
figure(3) 
plot(k,abs(fftshift(SIG)/fs))

xlabel('Frequency (Hz)')
ylabel('Magnitude')
title ('First received signal in frequency domain')

%First received signal in frequency domain (Rescaled to between 225Hz -250Hz)
figure(4) 
plot(k,abs(fftshift(SIG)/fs))
xlim([225 250])

xlabel('Frequency (Hz)')
ylabel('Magnitude')
title ('First received signal in frequency domain rescaled')


%% ------------3.4 Estimate the periodic noise--------------

 %Plot of the possible noise signal vs Signal 1 for 3 seconds
figure(5)
for i = 1:length(candidateT) 
    noise = estimateNoise(sig1,candidateT(i)); 
    noise = repmat(noise,[3 1]); 
    %plot Signal 1 
    subplot(5,2,i) 
    plot(t,sig1,'y') 
    hold on
    q = t(1:length(noise));  
    %plot noise signal 
    plot(q,noise,'r') 
    xlim([0,3]) 
    
    title(sprintf('Periodic Noise Signal %d Vs Time', i))
    ylabel('Magnitude')
    xlabel('Time (sec)')
end

% From figure 5 we can esimate signal 7 to be the noise signal as the shape
% of the signals closely match.


%Store noise profile 
Noisesig = estimateNoise(sig1,candidateT(7)); 

%Plot and compare noise signal to signal 1
figure(6)
plot(t(1:length(Noisesig)),sig1(1:length(Noisesig)),'y') 
grid on
hold on
plot(t(1:length(Noisesig)),Noisesig,'r','LineWidth',1.5)

title('Plot of selected Periodic Noise Signal (Signal #7)')
xlabel('Time (sec)')
ylabel('Magnitude')
legend ('Signal 1','Noise Signal 7','Location','southoutside')


%From figure 6 we determine that the period is bewteen 1.5 and 1.6 seconds
% Used T = candidateT(7)*Ts and got 1.531 
T = 1.531; 


%% -------------3.5 Model the periodic noise------------------

%Transpose Noisesig
Trans_Noisesig = Noisesig';

% Using Trigonometric Fourier Series 
f0 = 1/T;
a0 = (1/T).*sum(Trans_Noisesig.*Ts); 
h = t(1:candidateT(7));
FT_Noiesig = a0;

 
%For first 6 harmonics 
Harm = 10; %Chaanged to 10 harmonics for a much clearer image (Refer section 3.8)
for n = 1:Harm
    an = (2/T).*sum(Trans_Noisesig.*cos(2.*pi.*f0.*n.*h))*Ts;
    bn = (2/T).*sum(Trans_Noisesig.*sin(2.*pi.*f0.*n.*h))*Ts;
    FT_Noiesig = FT_Noiesig + an.*cos(2.*pi.*f0.*n.*h) + bn.*sin(2.*pi.*f0.*n.*h);
end

%Plot of Periodic noise and Modelled noise
figure(7) 
plot(t(1:length(Noisesig)),Noisesig,'b');
hold on
plot(t(1:length(FT_Noiesig)),FT_Noiesig,'r');
grid on

title('Periodic noise and Modelled noise')
xlabel('Time (seconds)')
ylabel('Magnitude')
legend ('Noise Signal', 'Modelled Noise','Location','southoutside')


%% ---------------------------3.6 Bias-----------------

a0 = 0;

%% -----------------3.7 Generate the approximation-----------------

Noisesig_fs = a0;


 %For 6 harmonics
for n = 1:Harm
    an = (2/T).*sum(Trans_Noisesig.*cos(2.*pi.*f0.*n.*h))*Ts;
    bn = (2/T).*sum(Trans_Noisesig.*sin(2.*pi.*f0.*n.*h))*Ts;
    Noisesig_fs = Noisesig_fs + an.*cos(2.*pi.*f0.*n.*h) + bn.*sin(2.*pi.*f0.*n.*h);
end

%Calculate number of times the signal was repeated and creates a vector
Repeat = ceil(length(t)/candidateT(7));
Noisesig_fs = repmat(Noisesig_fs, [1,Repeat]);
Noisesig_fs = Noisesig_fs(1:length(t));



%% ------------3.8 Compare the approximation----------------- 
% Plot received Signal 1 and New Modelled noise
figure(8) 
plot(t,sig1,'b');
hold on
plot(t,Noisesig_fs,'r');
grid on
xlim([0 3])

title(sprintf('Received Signal 1 and New Modelled noise with %d Harmonics',Harm))
xlabel('Time (seconds)')
ylabel('Magnitude')


%% ------------3.9 De-Noise---------------------------

%Removing periodic noise 
iml = sig;
im1 = sig1 - Noisesig_fs;

% Site 1 without periodic noise 
figure(9)
title('Site 1')
imshow(reshape(im1(1,:),480,640));



%% ----------------3.10 Remove the bandlimited random noise-------------

%Perform Fourier Transform to change signal to frequency domain.

IM1 = fftshift(fft(im1));
    
%remove the bandwidth random noise we calculated from section 3.3 
for i = 1:length(k)
   if ((k(i)>=200 && k(i)<=300)||(k(i)<=-200 && k(i)>=-300)) 
     
      IM1(i)= 0; 
   end
end

%Change final outcome to time domain
im2 = ifft(ifftshift(IM1));
    
%Display landing site 1 
figure (10)
title ("Image from signal 1")
imshow(reshape(im2,480,640));


%% -------------3.11 Choose a site---------------------

%Perform 3.9 and 3.10 for all four signals
for i = 1:4
    im1(i,:) = sig(i,:) - Noisesig_fs;
    IM1(i,:) = fftshift(fft(im1(i,:)));
end

%Apply filter for all four singals 
for i = 1:4
    SIGNAL1 = IM1(i,:);
    for x = 1:length(k)
        if ((k(x)>=200 && k(x)<=300) || (k(x)<=-200 && k(x)>=-300)) 
   
            SIGNAL1(x)= 0;
        end
    end
   
    IM1(i,:) = SIGNAL1;
    im2(i,:) = ifft(ifftshift(IM1(i,:)));
end

%Display final result 
for i = 1:4
    figure(11) 
    subplot(2,2,i)
    
    imshow(reshape(im2(i,:),480,640))
    title(sprintf('Site Number %d',i))
end


%% ----------------------3.12 Resolution----------------------------

% Display all images 
for i = 1:4
    figure(11+i)
    
    imshow(reshape(im2(i,:),480,640))
    title(sprintf('Site Number %d',i))
end




%% -------------------------------------------------------------------




