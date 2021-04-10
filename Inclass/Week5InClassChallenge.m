% Week 5 Code challenge
clear;
close all;
clc;
set(0,'DefaultFigureWindowStyle','docked');
% Note some common variable naming standards:
% Start with lower case letter for time domain signal
% i.e x, y
% Start with Upper case letter for frequency domain signal
% i.e. X, Y

%% Part A
% --- Q1 ---
% Define commonly used variables and magic numbers
samples = 2000;
T = 10;

% Time vector goes from 0 to 10 seconds
tA = linspace(0, T, samples+1);
tA(end) = [];

% --- Q2 ---
% Define commonly used variables and magic numbers
f = [5 14 17 40 63];
n = 1:5;

% Calculate the composite waveform
phase = repmat(n'*pi/4, [1 length(tA)]);
x = cos(2*pi*f'*tA + phase);
x = sum(x);

% --- Q3 --- 
% Find sampling frequency to be able to determine the frequency vector
fsA = 1/(tA(2) - tA(1));

% Calculate frequency vector
fA = linspace(-fsA/2, fsA/2, length(tA) + 1);
fA(end) = [];

% --- Q4 --- 
% Transform to frequency domain
X = fft(x);

% NOTE: Careful to note that if 'x' is a 1 x samples we can fft this straight 
% up as shown above but if 'a' is a matrix of signals, say 
% num_signals x samples, we must fft each row individually and NOT do fft(a)

% --- Q5 --- 
% NOTE:
% Divide by scaling factor when viewing the magnitude spectrum as MATLAB does
% scaling when using fft(). 
% Divide by length of signal if the signal is simple and made up of constant
% pure sinosoids
% (i.e. cos(...) + sin(...)
% Divide by fs if otherwise
% (i.e. a speech signal)

% Plot the signal in the time and frequency domains
figure
subplot(3,1,1)
plot(tA, x)
title("Time domain")

subplot(3, 1, 2)
plot(fA, fftshift(abs(X))/fsA)
title("Magnitude spectrum")

subplot(3, 1, 3)
plot(fA, fftshift(angle(X)))
title("Phase spectrum")

% Lets have a closer look at the time domain as it's hard to see what's
% happening


%% Part B

% --- Q1 ---
% Define commonly used variables and magic numbers

% Time vector goes from -50 to +50 seconds
tB = linspace(-50, 50, samples + 1);
tB(end) = [];

% --- Q2 --- 
% Define commonly used variables and magic numbers
fc = 5;
% Generate the signal and its shifted version
y = sinc(tB);
yshift = y.*cos(2*pi*fc*tB);

% --- Q3 --- 
% Find sampling frequency to be able to determine the frequency vector
fsB = 1/(tB(2) - tB(1));

% Calculate the frequency vector 
fB = linspace(-fsB/2, fsB/2, length(tB)+1);
fB(end) = [];

% --- Q4 --- 
% Transform to frequency domain
Y = fft(y);
YShift = fft(yShift);

% --- Q5 --- 
% Remember from before:
% Divide by scaling factor when viewing the magnitude spectrum as MATLAB does
% scaling when using fft(). 
% Divide by length of signal if the signal is simple and made up of constant pure sinosoids
% (i.e. cos(...) + sin(...)
% Divide by fs if otherwise
% (i.e. a speech signal)
figure

subplot(2, 2, 1)
plot(tB, y)
title("y(t)")

subplot(2, 2, 2)
plot(fB, fftshift(abs(Y))/fsB)
title("Y(f)")

subplot(2, 2, 3)
plot(tB, yShift)
title("yShift(t)")

subplot(2, 2, 4)
plot(fB, fftShift(abs(YShift))/fsB)
title("YShift(f)")