function [a0, an, bn, s_approx, T] = trigFS4(s_hinf, t, N)

    %N is an input
    %s_hinf is the function
    
%     nHarm = N;
%     n_trig = 1:nHarm;
%     f0 = 1/T;
%     Ts = t(2) - t(1);
%     T = Ts*length(t);
%     a0 = 1/T * sum(s_hinf) * Ts;
%     an = 2/T * s_hinf * cos(2*pi*f0*n_trig'*t).' * Ts;
%     bn = 2/T * s_hinf * sin(2*pi*f0*n_trig'*t).' * Ts;
%     s_approx = a0 + an * cos(2*pi*n_trig*t/T) + bn * sin(2*pi*n_trig*t/Tt);
    
%     Tstamp = t(2) - t(1); %Getting the time stamp
%     T = 1; %Period is one
%     n = (1:N)'; 
%     %T = Ts*length(t);
%     w0nt = 2*pi*n'*t/T;
%     a0 = 1/T * sum(s_hinf) * Tstamp;
%     an = 2/T * (s_hinf * cos(w0nt).') * Tstamp;
%     bn = 2/T * (s_hinf * sin(w0nt).') * Tstamp; 
%     s_approx = a0 + an * cos(w0nt) + bn * sin(w0nt);  
%     Tstamp = t(2) - t(1); 
%     T = t(end) - t(1) + Tstamp;
%     t = linspace(0,T,100+1); t(end) = [];
%     n = (1:N)'; 
%     
%     2*pi*f0*n_trig' *t
%     
%     w0nt = 2 * pi * n * t/T;
%     a0 = 1/T * sum(s_hinf) * Tstamp;
%     an = 2/T * (s_hinf * cos(w0nt).') * Tstamp;
%     bn = 2/T * (s_hinf * sin(w0nt).') * Tstamp; 
%     s_approx = a0 + an * cos(w0nt) + bn * sin(w0nt);  
%% Code i can submit
    
%    Test code
     t = linspace(0,1,100+1); t(end) = [];
     s_hinf = exp(t);
     N = 15;
     
    Tstamp = t(2) - t(1); 
    T = t(end) - t(1) + Tstamp;
    n = (1:N)';
    subfunc = 2 * pi * n * t * (1/T);
    
    a0 = 1/T * sum(s_hinf) * Tstamp;
    an = 2/T * s_hinf*cos(subfunc).'*Tstamp;
    bn = 2/T * s_hinf*sin(subfunc).'*Tstamp;
    
    s_approx = a0 + an * cos(subfunc) + bn * sin(subfunc);
    figure();
    plot(t, s_approx);
    
%% darcy code
%     Ts = t(2) - t(1);
%     T = t(end) - t(1) + Ts;
%     f0 = 1/T;
%     n_trig = 1:N;
%     
%     a0 = 1/T * sum(s_hinf) * Ts;
%     an = 2/T * s_hinf*cos(2*pi*f0*n_trig' *t).'*Ts;
%     bn = 2/T * s_hinf*sin(2*pi*f0*n_trig' *t).'*Ts;
%     
%     s_approx = a0 + an * cos(2*pi*f0*n_trig' *t) + bn * sin(2*pi*f0*n_trig' *t);
%     
%     figure();
%     plot(t, s_approx);
end