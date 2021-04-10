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
hold on
plot(t,noiseSound)

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

hold on
%Time period tPeriod
T = 5;
tPeriod = linspace(0,T, 100+1);tPeriod(end) = [];

%Plot the first noise which is the piecewise function i.e s3_t
s3_t = (A*tPeriod)/4;
s3_t( tPeriod>= 2.5 & tPeriod < 5) = (-A*tPeriod (tPeriod>=2.5 & tPeriod<5)+ 5*A)/4 ;
additive_noise_first = s3_t;
plot(tPeriod,additive_noise_first, 'r', 'LineWidth', 2);

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


    %Test code
    N = 5;
    t = linspace(0,1,100+1); t(end) = [];
    
    %Defining 
    Tstamp = t(2) - t(1); 
    T = t(end) - t(1) + Tstamp;
    n = (-N:N)';
    subfunc = 2 * pi * n * t * (1/T);
    
    %Coeficients
    c0 = 1/T * sum(additive_noise_second) * Tstamp;
    cn = 1/T * additive_noise_second * exp(1j*-subfunc).'*Tstamp;
    
    s_approx = cn * exp(1j*subfunc);
    figure()
    plot(tPeriod,s_approx)
    title('Graphs of evaluated coeficients using trignometric series')
    xlabel('Time period given by tPeriod')
    ylabel('Magnitude')



%% Question 5
%Hand calcs and then calculate for harmonics  | Use the trignometric series
%first period

approxq4 = zeros(1, 100);
a0 = 25/8;
an = -(25*(2*sin((pi*n)/2)^2 - n*pi*sin(pi*n)))/(8*n^2*pi^2) - (250*sin((pi*n)/2)^2 - 250*sin(pi*n)^2 + 125*n*pi*sin(pi*n))/(40*n^2*pi^2);
bn = (5*(5*sin(pi*n) - 5*n*pi*cos(pi*n)))/(8*n^2*pi^2) + (25*(sin(pi*n) - 2*cos(pi*n)*sin(pi*n) + n*pi*cos(pi*n)))/(8*n^2*pi^2);


for n = 1:5
    approxq4 = approxq4 + - (25*(2*sin((pi*n)/2)^2 - n*pi*sin(pi*n)))/(8*n^2*pi^2) - (250*sin((pi*n)/2)^2 - 250*sin(pi*n)^2 + 125*n*pi*sin(pi*n))/(40*n^2*pi^2)*cos(2 * pi * n .* t * (1/T)) + (5*(5*sin(pi*n) - 5*n*pi*cos(pi*n)))/(8*n^2*pi^2) + (25*(sin(pi*n) - 2*cos(pi*n)*sin(pi*n) + n*pi*cos(pi*n)))/(8*n^2*pi^2)*sin(2 * pi * n .* t * (1/T));
end

approxq4 = approxq4 + a0;

%% Question 6
%For time tPeriod;
    Tq6 = 5;
    N = 5;
    
    tPeriodq6 = linspace(0,Tq6, length(noiseSound)/4+1);tPeriodq6(end) = [];
    Tstamp_q6 = tPeriodq6(2) - tPeriodq6(1); 
    %Tq6 = tPeriodq6(end) - tPeriodq6(1) + Tstamp_q6;
    n = (-N:N)';
    subfunc_q6 = 2 * pi * n * tPeriodq6 .* (1/Tq6);  
    s3_t_q6 = (A*tPeriodq6)/4;
    s3_t_q6( tPeriodq6>= 2.5 & tPeriodq6 < 5) = (-A*tPeriodq6 (tPeriodq6>=2.5 & tPeriodq6<5)+ 5*A)/4 ;
    additive_noise_first_q6 = s3_t_q6;
    
    %Coeficients
    c0_q6 = 1/Tq6 * sum(additive_noise_first_q6) * Tstamp_q6;
    cn_q6 = 1/Tq6 .* additive_noise_first_q6 * exp(1j*-subfunc_q6).'*Tstamp_q6;
    FS1st = cn_q6' .* exp(1j*subfunc_q6);
    FS1st = sum(FS1st);

%s3_tq6 = (A*tPeriodq6)/4;
%s3_tq6( tPeriodq6>= 2.5 & tPeriodq6 < 5) = (-A*tPeriodq6 (tPeriodq6>=2.5 & tPeriodq6<5)+ 5*A)/4 ;
%FS1st = s3_tq6; 
plot(tPeriodq6,abs(FS1st))

%FS2nd = 

%% Question 7
% Get the full noise | use FSTotal
%FSTotal = FS1st + FS2nd;

%% Question 8

 dnSnd = FStotal' - noiseSound;
%% Question 9

%Mess around with n and see what happens 
%% Question 10

 %sound(noiseSound, 44100)
