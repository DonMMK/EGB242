function [c0, cn, s_approx, T] = complexFS(s_hinf, t, N)

%T =    
%c0 = 
%cn =
%s_approx =


% 
% clc; 
% T = 1;
% N = 1;
% t = linspace(0,T,100+1); t(end) = [];
% s_hinf = exp(t);
% Ts = t(2) - t(1);
% T = t(end) - t(1) + Ts;
% f0 = 1/T;
% n = (1:N).'; 
% w0nt = 2*pi*n*t/T;
% a0 = 1/T * sum(s_hinf)*Ts; 
% an = 2/T * (s_hinf * cos(w0nt).') * Ts; 
% bn = 2/T * (s_hinf * sin(w0nt).') * Ts;  
% c0 = a0;    
% cn_pos = (1/2) * (an - 1j * bn); 
% cn_neg = (1/2) * (an +1j * bn);   
% cn = [fliplr(cn_neg), c0, cn_pos];  
% n_new = (-N:N).';
% s_approx= cn * exp(1j * 2 * pi * f0 * n_new * t);   
% figure();
% plot(t, s_approx);
% grid on 
%% Code for me to use
    %Test code
    N = 1;
    t = linspace(0,1,100+1); t(end) = [];
    s_hinf = exp(t);
    
    %Defining stuff
    Tstamp = t(2) - t(1); 
    T = t(end) - t(1) + Tstamp;
    n = (-N:N)';
    subfunc = 2 * pi * n * t * (1/T);
    
    %Coeficients
    c0 = 1/T * sum(s_hinf) * Tstamp;
    cn = 1/T * s_hinf * exp(1j*-subfunc).'*Tstamp;
    
    s_approx = cn * exp(1j*subfunc);
    
    figure();
    plot(t, s_approx);
    
%% DARCy code

%     Ts = t(2) - t(1);
%     T = t(end) - t(1) + Ts;
%     f0 = 1/T;
%     Ts = t(2) - t(1);
%     n_comp = -N:N;
%     
%     c0 = 1/T * sum(s_hinf) * Ts;
%     cn = 1/T * s_hinf * exp(1j*-2*pi*f0*n_comp'*t).'*Ts;
    

end