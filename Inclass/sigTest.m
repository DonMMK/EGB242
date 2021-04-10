function [Type] = sigTest(s, Samp)
% EGB242: Week 3 In Class Challenge
% By: Thomas Nugent
%
% This function determines if an input function is even, odd or neither
% and returns its type as a string.
% [Type] = sigTes(s, Samp)
%
% Will only work for signals whose first period is centred on zero.

% Set a threshhold for acceptable error
thresh = 1e-12;

% Test to see if we can test a whole period
if length(s) < Samp
    error('ERROR: s1 is less than one period.');
end

% Remove first sample so that s is centred about zero
s = s(2:Samp);

% Test for even-ness
% s(t) = s(-t) => s(t) - s(-t) = 0
if (s - fliplr(s)) < thresh
    Type = 'even';
% Test for odd-ness
% s(t) = -s(-t) => s(t) + s(-t) = 0
elseif (s + fliplr(s)) < thresh
    Type = 'odd';
else
    Type = 'neither';
end

end