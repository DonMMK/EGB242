%% Assignment 1 - Part A: GenerateDataAssignment1A.m
%  Clearing workspace
clear; close all; clc;
%  This file generates the assignment data needed to complete Assignment 1 - Part A. 
%  Data is generated based on your student number.
%  This script only needs to be executed ONCE.
%
%  Make sure that the following files are in the same directory.
%  The directory path should not have spaces nor special characters.
%  1) data
%  2) initialise1A.p
%  3) GenerateDataAssignment1A.m
%  4) mission.m
%
%  Enter your student number into line 22.
%  Student numbers have the format "n01234567".
%  Omit the leading 'n' and leading '0', and enter it as 1234567.
%  It should only be a 7 or 8 digit number.
%  8 digit numbers will be truncated.

% Student numbers
sid = 10496262;

% Check student id length
if (length(num2str(sid)) == 7)||(length(num2str(sid)) == 8)
    %  Generate test signals and initialise the assignment variables.
    % Generate the noisy speech signal and 
    % the scale factor and DC shift of the noise signal before affecting the original speech signal
    [s1t, s2t, noiseSound] = initialise1A(sid); % do not modify this line.
    fprintf('\nScript ran successfully.\n')
else
    fprintf('Invalid student number.\nPlease enter a valid student number.\n');
end

%  GenerateDataAssignment1A.m file will load and work with the Data1A file directly.
%  You may close this script once Data1A has been saved.
%  This script does not need to be executed again.
