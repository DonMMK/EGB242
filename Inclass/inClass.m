clear; clc; close all;

%% Part A 
T=10; % Time for 1 period 
samples = 100; % Number of samples per period

t1 = linspace(-T/2, T/2, samples + 1); t1(end) = [];

y1 = sin(pi * t1)./(pi * t1);
%y1( t >= 0 && t<= T/4) = 0.2; Creating a another function 
y1(t1 == 0 ) = 1;

% y1 = sinc(t1);

y2 = (t1/10).^3;

figure
subplot(2,1,1)
plot(t1, y1) 

subplot(2,1,2)
plot(t1,y2)

% Repeat the period 10 times
s1 = repmat(y1, [1 10]);
s2 = repmat(y2, [1 10]);

t2 = linspace(-T/2, 10*T - T/2, samples * 10 + 1); t2(end) = [];

figure
subplot(2,1,1)
plot(t2, s1)
subplot(2,1,2)
plot(t2, s2)

%% Part B

function [a0,an,bn] = FS(s,Harm, f0)
    T = 1/f0;
    samples = length(s);
    
    t = linspace(0, T, samples + 1); t(end) = [];
    
    Ts = t2(2) - t(1); % Time Step
    Fs = 1/Ts;         % Sampling rate 
    
    a0 = 1/T * sum(s) * Ts;
    
    an = 2/T * s * cos(2*pi*f0*n*t).' * Ts;
    bn = 2/T * s * sin(2*pi*f0*n*t).' * Ts;
    
end
