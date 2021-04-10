function [a0, an, bn, s1t_approx, sid] = trigFS(sid)
    syms t, syms n, pi = sym('pi');
      sid = 10496262;
    % As a symbolic expression, save the expressions 
    % for the trigonometric coefficients below:
    a0 = 9/32;
    an = (1/(2*pi*n))*((3/4)*sin((6*pi*n)/4) + (cos((6*pi*n)/4))/(2*pi*n) - 1/(2*pi*n));
    bn = (-3/4)*((cos((3*pi*n)/2))) + (sin((3*pi*n)/2))/(2*(pi^2)*(n^2));
    % Find the trigonometric FS, evaluate the FS with
    % 2 harmonics and save the expression below:
    a1 = (1/(2*pi*1))*((3/4)*sin((6*pi*1)/4) + (cos((6*pi*1)/4))/(2*pi*1) - 1/(2*pi*1));
    a2 = (1/(2*pi*2))*((3/4)*sin((6*pi*2)/4) + (cos((6*pi*2)/4))/(2*pi*2) - 1/(2*pi*2));
    b1 = (-3/4)*((cos((3*pi*1)/2))) + (sin((3*pi*1)/2))/(2*(pi^2)*(1^2));
    b2 = (-3/4)*((cos((3*pi*2)/2))) + (sin((3*pi*2)/2))/(2*(pi^2)*(2^2));
    
    s1t_approx = a0 + a1*cos(2*pi*1*t) + a2*cos(2*pi*2*t) + b1*sin(2*pi*1*t) + b2* sin(2*pi*2*t) ;
    
end

