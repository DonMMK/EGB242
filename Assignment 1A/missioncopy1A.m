%% Assignment 1 Part A - Section A2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%.%
%  Do not change before line 28.
%  If you have not generated Data1A from GenerateDataAssignment1A.m,
%  do that now.

%  Clearing and preparing the workspace
clear; clc; close all;

%  Load assignment data.
load Data1A;  

% VARIABLES:
% t - Time vector
% T - Period
% additive noise - Your noise waveform
% a0, an, bn - Trig Fourier series variables
% OR
% c0, cn - Complex Fourier series variables
% FS1 - Fourier series approximation
% dnSnd - De-noised resulting wave

%==================================================================
% Refer to the assignment sheet for details on variable naming.
% Names of the variables are important,
% e.g. 'a1' is considered a different variable to 'A1'.
%====Enter your code below this line================================
%% Question 1
%Coeficients
A = 5; B = 3;

%Plot the noiseSound
t = linspace(0,20,882000+1);t(end) = [];
figure();
plot(t,noiseSound)
hold on
%Plot sound one
T1 = 5;
t1 = linspace(0,T1, 100+1);t1(end) = []; 
t1 = repmat(t1, [1 4]);
s1_t = t1/A;
t1 = linspace(0, 4*T1, 400+1); t1(end) = [];
plot(t1,s1_t, 'g', 'LineWidth', 2)

%Plot sound two
T2 = 5;
t2 = linspace(0,T2, 100+1);t2(end) = []; 
t2 = repmat(t2, [1 4]);
s2_t = exp(-(t2-B)/4);
t2 = linspace(0, 4*T2, 400+1); t2(end) = [];
plot(t2,s2_t, 'b', 'LineWidth', 2)

%Plot sound three
t0_3 = 0;
T3 = 5;
t3 = linspace(t0_3 , t0_3 +T3, 100 + 1); t3(end) = [];
t3 = repmat(t3, [1 4]);
s3_t = (A*t3)/4;
s3_t( t3>= 2.5 & t3 < 5) = (-A*t3 (t3>=2.5 & t3<5)+ 5*A)/4 ;
t3 = linspace(t0_3 , t0_3 + 4*T3, 400 + 1); t3(end) = [];
plot(t3,s3_t, 'r', 'LineWidth', 2);
title('Plot of Received signals and the possible signals vs time');
xlabel('Time');
ylabel('Magnitude');

hold off
%% Question 2

%Time period tPeriod
T = 5;
tPeriod = linspace(0,T, 100+1);tPeriod(end) = [];
figure();
%Plot the first noise which is the piecewise function i.e s3_t
s3_t = (A*tPeriod)/4;
s3_t( tPeriod>= 2.5 & tPeriod < 5) = (-A*tPeriod (tPeriod>=2.5 & tPeriod<5)+ 5*A)/4 ;
additive_noise_first = s3_t;
plot(tPeriod,additive_noise_first, 'r', 'LineWidth', 2);


hold on
%Plot the second noise which is the linear function i.e s1_t

additive_noise_second = tPeriod/A;
plot(tPeriod, additive_noise_second, 'g', 'LineWidth', 2)

title('Graphs of interfering waveforms')
xlabel('Time period given by tPeriod')
ylabel('Magnitude')
hold off
%% Question 3

%To be done using the deifitions
%It's odd

%% Question 4
%Evaluate additive noise second using matlab linear func

    %Test code
    N = 5;
    t = linspace(0,1,100+1); t(end) = [];
    
    %Defining 
    Tstamp = t(2) - t(1); 
    %T = t(end) - t(1) + Tstamp;
    T = 5;
    n = (-N:N)';
    subfunc = 2 * pi * n * t * (1/T);
    
    %Coeficients
    c0 = 1/T * sum(additive_noise_second) * Tstamp;
    cn = 1/T * additive_noise_second * exp(1j*-subfunc).'*Tstamp;
    
    FS2nd = cn * exp(1j*subfunc);
    figure()
    plot(tPeriod,FS2nd)
    title('Approximate for noise second')
    xlabel('Time period given by tPeriod')
    ylabel('Magnitude')



%% Question 5
%Hand calcs and then calculate for harmonics  | Use the trignometric series
%first period
% T= 5;
% approxq4 = zeros(1, 100);
% a0 = 25/16;
% an = ((-0.316629 + 0.316629*cos(3.14159.*n) + 0.994718.*n*sin(3.14159.*n))/(n^2))  + ((0.0633257*cos(3.14159.*n) - 0.0633257*cos(6.28319.*n) - 0.198944.*n*sin(3.14159.*n) - (5.551115123125783^-17).*n*sin(6.28319*n))/(n.^2))  ;
% bn = ((-0.994718.*n*cos(3.14159.*n) + 0.316629*sin(3.14159.*n))/(n.^2) )  + ((0.0633257*cos(3.14159.*n) - 0.0633257*cos(6.28319.*n) - 0.198944.*n*sin(3.14159.*n) - (5.551115123125783^-17).*n*sin(6.28319.*n))/(n.^2));
% 
% %To get for 5 harmonics
% for n = 1:5
%     approxq4 = approxq4 + an*cos(2 * pi .* n .* t * (1/T)) + bn*sin(2 * pi * n .* t * (1/T));
% end
% 
% approxq4 = approxq4 + a0;

    %Test code
    N = 5;
    t = linspace(0,1,100+1); t(end) = [];
    
    %Defining 
    Tstamp = t(2) - t(1); 
    %T = t(end) - t(1) + Tstamp;
    T = 5;
    n = (-N:N)'; n(n == 0) = [];
    subfunc = 2 * pi * n * t * (1/T);
    
    %Coeficients
    c0_q5 = 1/T * sum(additive_noise_first) * Tstamp;
    cn_q5 = 1/T * additive_noise_first * exp(1j*-subfunc).'*Tstamp;
    
    FS1st = cn_q5 * exp(1j*subfunc);
    figure();
    plot(tPeriod,FS1st)
    title('Approximate for noise first')
    xlabel('Time period given by tPeriod')
    ylabel('Magnitude')


%% Question 6

 %FS1st = cn .* exp(1j*subfunc);
 
 %FS2nd = cn_q5 .* exp(1j*subfunc);

%% Question 7
% Get the full noise | use FSTotal
 FSTotal = FS1st + FS2nd;

%% Question 8

dnSnd = FStotal - noiseSound;
%% Question 9
%Mess around with n and see what happens 
%% Question 10
 %sound(noiseSound, 44100)
