%% Question 4

A = 5;
%TQ4 = 5;
f0 = 1/TQ4;
N = 5;
n = (1:N).';
Tstamp = tPeriod(2) - tPeriod(1);
TQ4 = tPeriod(end) - tPeriod(1) + Tstamp;
tPeriod = linspace(0,TQ4,(44100*5) +1);tPeriod(end) = [];
additive_noise_second = tPeriod/A;

% Tstamp = t(2) - t(1); 
% %T = t(end) - t(1) + Tstamp;
% %n = (1:N)';
% subfunc = 2 * pi * N' .* t * (1/T);
%     
% a0 = 1/T .* sum(additive_noise_second) .* Tstamp;
% an = 2/T .* additive_noise_second.*cos(subfunc).'*Tstamp;
% bn = 2/T .* additive_noise_second.*sin(subfunc).'*Tstamp;
% 
% sapprox = a0 + an.*cos(2*pi*N'.*tPeriod/T) + bn.*sin(2*pi*N'.*tPeriod/T);

    subfunc = 2 * pi * n * f0 * tPeriod;
    
    a0 = 1/TQ4 * sum(additive_noise_second) * Tstamp;
    an = 2/TQ4 * additive_noise_second*cos(subfunc).'*Tstamp;
    bn = 2/TQ4 * additive_noise_second*sin(subfunc).'*Tstamp;
    
    s_approx = a0 + an * cos(subfunc) + bn * sin(subfunc);
    
    figure();
    
    plot(tPeriod, s_approx);