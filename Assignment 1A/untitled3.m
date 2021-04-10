clc

syms t n

T = 5; %T is the period

f0 = 1/T; % f0 is the fundamental frequency

a0 = (1/T)*(int((5*t)/4,t,0,2.5) + int((5*t)/4,t,2.5,5))

an_1 = (2/T)*(int((5*t/4)*cos(2*pi*n*f0*t),t,0,2.5));

an_2 = (2/T)*(int(((-5*t + 25)/4)*cos(2*pi*n*f0*t),t,2.5,5));

an_FS1st = an_1 + an_2
 
bn_1 = (2/T)*(int((5*t/4)*sin(2*pi*n*f0*t),t,0,2.5));

bn_2 = (2/T)*(int(((-5*t + 25)/4)*sin(2*pi*n*f0*t),t,2.5,5));

bn_FS1st = bn_1 + bn_2