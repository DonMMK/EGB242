t = linspace(-1, 1, 2^12 + 1); t(end) = [];
s2t = 8*rectangularPulse(-5/2, 5/2, 1)*cos(60*pi*t);

% S2f = fft(s2t);

S2f = 20*(sinc(5*t - 30) + sinc(5*t + 30));

figure();
% plot(t, s2t);
% hold on
plot(t, S2f);