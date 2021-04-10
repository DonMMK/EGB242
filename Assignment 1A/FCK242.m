clc

syms t n

T =2; %T is the period

f0 = 1/T; % f0 is the fundamental frequency

%n=5;

a0 = (1/T)*(int(t,t,0,3/4) + int(0,t,3/4,1))

%an = (2/T)*(int(t*cos(2*pi*n*f0*t),t,0,3/4) + (int(0*cos(2*pi*n*f0*t),t,3/4,1)))
 
%bn = (2/T)*(int(t*sin(2*pi*n*f0*t),t,0,3/4) + (int(0*sin(2*pi*n*f0*t),t,3/4,1)))

 c0 = (1/T)*( (int((exp(3*t) ),t,-0.5,0)) + (int((exp(-3*t) -2 ),t,0,0.5)) )
 
  cn = (1/T)*( (int((exp(3*t))*(exp(-2*pi*i*n*t)),t,-0.5,0)) + (int((exp(-3*t) -2 )*(exp(-2*pi*i*n*t)),t,0,0.5)) )

% cn = (1/T)*int((exp(3*t) - 2 )*(exp(-1j*2*pi*n*t)),t,0,1)
 
%cn = (1/(3 - 2*i*pi*n*1)) - (  exp(1.5 + i*pi*n*1)  / (3 - 2*i*pi*n*1)) - (  exp(-1.5 - i*pi*n*1) /  (3 + 2*i*pi*n*1)  ) + ( 1/(3 + 2*i*pi*n*1) ) + ((exp(-i*pi*n*1))/(i*pi*n*1)) - (1/ (i*pi*n*1))
% s2t = (cn.*exp(1j*2*pi*n*t))

 %sid = 10496262;
%[c0, cn, s2t_approx, sid] = compFS(sid);

% This line will output the two functions generated through generateDataAssignment1A.m
%[s1t, s2t] = generateFunction(sid)

% This line will display a typeset version of your input solution
%displayResultQ2(c0, cn, s2t_approx)   