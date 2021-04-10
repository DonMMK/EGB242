%% Tutorial 10 In Class Challenge
clear all, close all, clc,

%% Question 1 - ilaplace and laplace

% 1. 
syms s;
invLap1 = ;
pretty(invLap1)

% 2. 
syms ;
invLap2 = ;
pretty(invLap2)

% 3. 
syms ;
cosFunc = ;
CosFunc = ;
cosFuncFinal = ;
pretty(cosFunc); pretty(CosFunc); pretty(cosFuncFinal);

%% Question 2 - Square Brackets and Numerical Arrays

% 1. 
help randn
test1 = ;

% 2.
test2 = ; % allocate memory
tstart = ; % use tic
for 


% 3. 

test3 = ;


% 4. 
help mean
timeAv2 = ;
timeAv3 = ;

ensAv2 = ;
ensAv3 = ;

figure(1), clf
subplot(2,1,1), plot(timeAv2), hold on, plot(timeAv3)
subplot(2,1,2), plot(ensAv2), hold on, plot(ensAv3)
% 
subplot(2,1,1), plot((1:25)/5,'r*')
subplot(2,1,2), plot(mean((1:25)/5)*ones(1,1e3),'r*')

%% Question 3 - Curly Braces and Cells

% Create a cell array with a numerical array in the first element, a string
% in the second and a single number in the third
cellArray = {};

% 1. Print the second element to the command window
disp();

% 2. Open the cell in the workspace and navigate to the third element
% What can we say about accessing the element and navigating to it? When
% you navigate MATLAB is just accessing whatever you click on

% 3. Create a cell with a ?for? loop over 1 to 5 and store 1*s, 2*s ... all the way
% up to 5*s in the cell. Investigate how it is stored in MATLAB.
syms s;
for i = 1:5
    
end

%% Question 4

% 1. Create a structure with a numerical array in the first element, a string in the
% second and a cell array in the third. Give them meaningful field names.
structure = struct();

% 2. Print the second element of the numerical array to the command window.
disp();

% 3. Create the transfer function displayed in the first part of this tute as a
% transfer function object in MATLAB.
tfunc = tf();

% 4. Replace the numerator with the equivalent numerical array meaning s^2 + 1
tfunc. = ;
% OR
tfunc. = ;
% OR
tfunc. = ;