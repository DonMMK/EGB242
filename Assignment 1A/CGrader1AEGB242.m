%% Chirans code

%% S1Q3
clc
close all

T = 1;
t = linspace(0, 1*T, 100+1); t(end) = [];
t = repmat(t, [1 5]);
s2_hinf = exp(1 - 2*t) + 6;
t = linspace(0, 5*T, 500+1); t(end) = [];

figure()
plot(t,s2_hinf)

%% S1Q4
T = 1;
N = 2;
t = linspace(0,T,100+1); t(end) = [];
s_hinf = exp(1-(2*t))+6;
n = (1:N)';

Ts = t(2) - t(1);

T = Ts*length(t);

w0nt = 2*pi*n*t/T;

a0 = 1/T * sum(s_hinf) * Ts;

an = 2/T * (s_hinf * cos(w0nt).') * Ts;

bn = 2/T * (s_hinf * sin(w0nt).') * Ts; 

s_approx = a0 + an * cos(w0nt) + bn * sin(w0nt);

figure();
plot(t,s_approx,'red--');
grid on

%% S1Q5
clc; 
T = 1;
N = 2;
t = linspace(0,T,100+1); t(end) = [];
Ts = t(2) - t(1);
T = Ts*length(t);
f0 = 1/T;
n = (1:N).'; 
w0nt = 2*pi*n*t/T;
s_hinf = exp(1-(2*t))+6;
a0 = 1/T * sum(s_hinf)*Ts; 
an = 2/T * (s_hinf * cos(w0nt).') * Ts; 
bn = 2/T * (s_hinf * sin(w0nt).') * Ts;  
c0 = a0;    
cn_pos = (1/2) * (an - 1j * bn); 
cn_neg = (1/2) * (an +1j * bn);   
cn = [fliplr(cn_neg), c0, cn_pos];  
n_new = (-N:N).';
s_approx= cn * exp(1j * 2 * pi * f0 * n_new * t);   
figure();
plot(t, s_approx);
grid on