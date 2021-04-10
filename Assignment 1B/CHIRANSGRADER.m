%% Question 1
function [sid,S1f] = fourierT(sid)
    % Do not alter line below:
    syms f
    % Enter your student ID as sid:
    sid = 10454012;
    % As a symbolic expression, save the expressions 
    % for S1f below:
    S1f = 7/(9 + 1j*2*pi*f)
end

%% Question 2
function [sid,S2f] = fourierT(sid)
    % Do not alter line below:
    syms f
    % Enter your student ID as sid:
    sid = 10454012;
    % As a symbolic expression, save the expressions 
    % for S2f below:
    S2f = 9/2 * ((1 / (6 + 1j*2*pi*(f - 25))) + (1 / (6 + 1j*2*pi*(f + 25))));
end

%% Question 3
function [sid,S3f] = fourierT(sid)
    % Do not alter line below:
    syms t f
    % Enter your student ID as sid:
    sid = 10454012;
    % As a symbolic expression, save the expressions 
    % for S3f below:
    S3f = 1/2*( (9/4*(1/(6 + 1j*2*pi*(f - 50)) + 1/(6 + 1j*2*pi*f)) + 7*sinc(7*(f - 25)))  +  (9/4*(1/(6 + 1j*2*pi*f) + 1/(6 + 1j*2*pi*(f + 50))) + 7*sinc(7*(f + 25))) );
end

%% Question 4
function [sid,f1,S1f] = fourierT(sid)
    % Enter your student ID as sid:
    sid = 10454012;
    % Define the frequency vector
    f1 = linspace(-50, 50, 10000 + 1); f1(end) =[];
    % evaluate S1f (do not use fft - use answer from previous questions)
    S1f = 7./(9 + 1j*2*pi*f1);
end

%% Question 5
function [sid,t2,s2t,Ts,fs,k2,S2k, mag_S2k, phase_S2k]  = q5(sid);
    % Enter your student ID as sid:
    sid= 10454012;
    % Generate Time Vector 
    t2 = linspace(-2, 2, 999 + 1); t2(end) = [];
    % Generate signal vector s2(t) in time domain
    s2t = 9.*exp(-6*t2).*heaviside(t2).*cos(50*pi*t2);
    % Determine the sampling period 
    Ts = t2(2) - t2(1);
    % Determine the Sampling Frequency
    fs= 1/Ts;
    % Generate frequency vector
    k2 = linspace(-fs/2, fs/2, 999 + 1); k2(end) =[];
    % Perform Discrete Fourier Transform of s2t
    S2k= fft(s2t);
    
    mag = abs(S2k);
    phase = angle(S2k);
    
    % store the magnitude plot as 'mag_S2k'
    figure(1)
    mag_S2k = plot (k2, fftshift(mag)/fs, 'blue');
    xlabel("Frequency");
    ylabel("Magnitude");
    title("Magnitude Spectrum of S2k");
    
    % store the phase plot as 'phase_S2k'
    figure(2)
    phase_S2k = plot (k2, fftshift(phase), 'blue');
    xlabel("Frequency");
    ylabel("Phase");
    title("Phase Spectrum of S2k");
end

%% Question 6
function [sid,t3,s3t,Ts,fs,k3,S3k, mag_S3k, phase_S3k] = q6(sid);
    % Enter your student ID as sid:
    sid= 10454012;
    % Generate Time Vector 
    t3 = linspace(-2, 2, 999 + 1); t3(end) = [];
    % Generate signal vector s3(t) in time domain
    s3t = (4.5 .* exp(-6*t3) .* heaviside(t3) .* cos(50*pi*t3) + rectangularPulse(-7/2,7/2,t3)) .* cos(50*pi*t3);
    % Determine the sampling period 
    Ts = t3(2) - t3(1);
    % Determine the Sampling Frequency
    fs= 1/Ts;
    % Generate frequency vector
    k3 = linspace(-fs/2, fs/2, 999 + 1); k3(end) =[];
    % Perform Discrete Fourier Transform of s3t
    S3k= fft(s3t);
    
    mag = abs(S3k);
    phase = angle(S3k);
    
    % store the magnitude plot as 'mag_S3k'
    figure(1)
    mag_S3k= plot (k3, fftshift(mag)/fs, 'blue');
    xlabel("Frequency");
    ylabel("Magnitude");
    title("Magnitude Spectrum of S3k");
    
    % store the phase plot as 'phase_S3k'
    figure(2)
    phase_S3k= plot (k3, fftshift(phase), 'blue');
    xlabel("Frequency");
    ylabel("Phase");
    title("Phase Spectrum of S3k");
end

%% END ??