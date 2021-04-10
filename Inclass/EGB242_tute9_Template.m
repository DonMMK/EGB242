% EGB242 Tutorial 9 Computer Lab
clear all, close all, clc

% Type 'doc tf'

% Set up numerator and denominator as per help file
num = [1 0 0];
den = [1 le3 le9;
% Create transfer function
tfunc = tf(num, den);

%% Question 1
% Poles and zeros
poles = pole(tf);
zeros = zero(tf);
% Display to command window
disp(poles), disp(zeros)

% Plot
figure(1), clf, pzmap(tfunc) % pole zero map

%% Question 2
figure(2),clf,  % bode plot


%% Question 3
ltiview(tfunc)
% LTIviewer provides a graphical interface to work with a function 

% Change line to red and dashed - 2 ways
% 1 - 
% 2 - Command line

ltiview(tfunc, '--r')
% Save to local disk
% 
% Other options

%% Question 4
% Characteristics of impulse response
3.16e4 * 2*pi

%% Question 5
% Set up new transfer function

