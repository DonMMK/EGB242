%% TEMP FUNCTION

%% Question 1
function [sid,S1f] = fourierT(sid)
    % Do not alter line below:
     syms f
    % Enter your student ID as sid:
    sid = 10496262;
    % As a symbolic expression, save the expressions 
    % for S1f below:
    
    f = linspace(-1 , 1, 2^12 + 1); f(end) = [];
    S1f = 2/(2 + 1j*pi*f);
    
    
    t = linspace(-1, 1, 2^12 + 1); t(end) = [];
    inv1 = ifft(S1f, 'symmetric');
    inv2 = ifft(S1f, 'nonsymmetric');
    ori = 4*exp(-4*t) * heaviside(1);
     
    figure()
    plot(f, S1f)
    
    figure()
    hold on
    plot(t, inv, 'r')
    
    plot(f, S1f, 'b');
    plot(t, inv1, 'r');
    plot(t, inv2, 'g');
    plot(t, ori, 'm');
    
end

%% Question 2
function [sid,S2f] = fourierT(sid)
    % Do not alter line below:
    syms f
    % Enter your student ID as sid:
    sid = 10496262;
    % As a symbolic expression, save the expressions 
    % for S2f below:
    S2f = 20*( sinc(5*(f - 30) ) + sinc(5*(f + 30) )  );
end

%% Question 3
function [sid,S3f] = fourierT(sid)
    % Do not alter line below:
    syms t f
    % Enter your student ID as sid:
    sid = 10496262;
    % As a symbolic expression, save the expressions 
    % for S3f below:
    S3f = 0.5*((2/(8 + 1j*2*pi*(f- 30)))+ 20*( sinc(5*(f - 60) ) + sinc(5*f )  ) +  (2/(8 + 1j*2*pi*(f+ 30))) + 20*( sinc(5*f ) + sinc(5*(f + 60) )  )     );
end

%% Question 4
function [sid,f1,S1f] = fourierT(sid)
    % Enter your student ID as sid:
    sid = 10496262;
    % Define the frequency vector
    f1 = linspace(-50, 50, 10000 + 1); f1(end) =[];
    % evaluate S1f (do not use fft - use answer from previous questions)
    S1f = 2./(2 + 1j*pi*f1);
end

%% Question 5
function [sid,t3,s3t,Ts,fs,k3,S3k, mag_S3k, phase_S3k] = q6(sid);
    % Enter your student ID as sid:
    sid= 10496262;
    % Generate Time Vector 
    t3 = linspace(-2,2, 999 + 1);
    t3(end)=[];
    % Generate signal vector s3(t) in time domain
    s3t = (( 2.*exp(-8*t3) .* heaviside(t3) + 8 .* rectangularPulse(-5/2 , 5/2 ,t3) .* cos(60*pi*t3)) .* cos(60*pi*t3)) ; % Check this function
    % Determine the sampling period 
    Ts = t3(2) - t3(1) ;
    % Determine the Sampling Frequency
    fs= 1/Ts;
    % Generate frequency vector
    k3 =linspace(-fs/2 , fs/2, 999 + 1);
    k3(end)=[];
    % Perform Discrete Fourier Transform of s3t
    S3k= fft(s3t) ;

    magnitude = abs(S3k);
    phase = angle(S3k);
    
    % storing the magnitude as 'mag_S3k'
    figure(1)
    mag_S3k= plot(k3, fftshift(magnitude)/fs , 'blue');
    xlabel('frequency')
    ylabel('magnitude')
    title('Plotting Magnitude');
    % storing the phase as 'phase_S3k'
    figure(2)
    phase_S3k= plot(k3, fftshift(phase),'blue') ;
    xlabel('frequency')
    ylabel('phase')
    title('Plotting Phase');
end