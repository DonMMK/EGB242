%% Tutorial 10 In Class Challenge
clear all, close all, clc,

%% Question 1 - ilaplace and laplace

% 1. 
syms s;
F2s = 10*(s+2)/(s*(s^2-4*s-5));
invLap1 = ilaplace(F2s);
pretty(invLap1)

% 2. 
syms s t;
invLap2 = ilaplace(1, s, t);
pretty(invLap2)

% 3. 
syms a t;
cosFunc = cos(a*t);
CosFunc = laplace(cosFunc);
cosFuncFinal = ilaplace(CosFunc);
pretty(cosFunc); pretty(CosFunc); pretty(cosFuncFinal);

%% Question 2 - Square Brackets and Numerical Arrays

% 1. 
help randn
test1 = randn(1, 1000);

% 2.
itter = 25;
test2 = ones(itter, 1000); % allocate memory
tstart = tic(); % use tic
for ii = 1:itter
    test2(ii,:) = ii/5 + 10/ii*randn(1,1000);
end
toc(tstart)


% 3. 
tstart2 = tic(); % use tic
test3 = (1:itter).'/5 + 10./(1:itter).'.*randn(itter,1000);
toc(tstart2)


% 4. 
help mean
timeAv2 = mean(test2);
timeAv3 = mean(test3);

ensAv2 = mean(test2);
ensAv3 = mean(test3, 1);

figure(1), clf
subplot(2,1,1), plot(timeAv2), hold on, plot(timeAv3)
subplot(2,1,2), plot(ensAv2), hold on, plot(ensAv3)
% 
subplot(2,1,1), plot((1:25)/5,'r*')
subplot(2,1,2), plot(mean((1:25)/5)*ones(1,1e3),'r*')

%% Question 3 - Curly Braces and Cells

% Create a cell array with a numerical array in the first element, a string
% in the second and a single number in the third
cellArray = {[1 2 3], 'string', 7};

% 1. Print the second element to the command window
disp(cellArray{2} (2));

% 2. Open the cell in the workspace and navigate to the third element
% What can we say about accessing the element and navigating to it? When
% you navigate MATLAB is just accessing whatever you click on

% 3. Create a cell with a ‘for’ loop over 1 to 5 and store 1*s, 2*s ... all the way
% up to 5*s in the cell. Investigate how it is stored in MATLAB.
syms s;
for i = 1:5
    cellArray2{i} = i*s;
end

%% Question 4

% 1. Create a structure with a numerical array in the first element, a string in the
% second and a cell array in the third. Give them meaningful field names.
structure = struct('array', [21 65 33], 'string', 'hello', 'cell', {cellArray});

% 2. Print the second element of the numerical array to the command window.
disp(structure.array(2));

% 3. Create the transfer function displayed in the first part of this tute as a
% transfer function object in MATLAB.
num = [1 0 0];
den = [1 10^3 10^9];
tfunc = tf(num, den);

% 4. Replace the numerator with the equivalent numerical array meaning s^2 + 1
tfunc.numerator = [1 0 1];
% OR
tfunc.numerator = {[1 0 1]};
% OR
tfunc.numerator{1}(3) = 1;