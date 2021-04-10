%Restting script
close all; clear all; clc

%% Take Home Computer Lab
%% Question 1
%My first Matlab script
T=1;
f= linspace(-0.5,0.5,100); %f = linspace(-0.5,0.5,samples+1); f = f(1:end-1);
y = (sin(pi*f*T))/(pi*f); 

figure
plot(f,y);
%title('frequency vs y ')
%xlabel('frequency')
%ylabel('y')


%% Question 2
%function [output1,output2,...] = FunctionName(input1,input2,...)
%check this again
%function [y] = myfunc(T)
%    samples = 100;
%    f = linspace(-0.5,0.5,samples+1);
%    f = f(1:end-1);
%    y = T*sin(pi*f*T)./(pi*f*T);
%    figure(1)
%    plot(f,y)
%end

%% Video Tutorial 
%% Question 5
A = 240;
f = 50;
T = 1/f;

%build time vector
t = linspace(0,5*T,2000); %(start,stop,number of points)

y = A*cos(2*pi*f*t);

figure(1)
plot(t,y)

y1 = 120*(exp(j*100*pi*t) + exp(-j*100*pi*t));
figure(2)
plot(t,y1)

%% Code Challenge
function [st, t] = gencos(f0, phi, R)
%UNTITLED Summary of this function goes here
%File needs to be called the same thing
%   Detailed explanation goes here
%Checking that R is an integer
if (round(R) ~= R)
    error('R is not an integer.');
end

%100 samples per the period
samples = 100;
T = 1/f0; %period of function

%Always use samples + 1 and remove the last element.
%create time vector t
t = linspace(0, T* R, samples* R + 1); 
%REMOVE LAST ELEMENT
t(end) = [];

st = cos(2*pi*f0*t + phi);

if R > 1
    randPoint() = round(rand() * samples);
    %if sum(st(1:samples) - st(samples + 1: samples*2)) < 1e-10
    if st(randPoint) - st(randPoint + samples) < 1e-10
        fprintf('Yes the function is periodic.\n')
    else
        fprintf('No the function is not periodic\n')
    end
else
    fprintf('Cannot determine if function is periodic because R is equal to 1')
end

figure
plot(t,st)
title(['Plot of the function cos(2*pi)',num2str(f0), 't + ' num2str(phi) + ')'])
xlabel('time')
ylabel('amplitude')


end
