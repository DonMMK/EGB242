function [c0, cn, s2t_approx, sid] = compFS(sid)
    syms n t; e = sym(exp(sym(1))); pi = sym('pi');
    % Enter your student ID as sid:
    sid = 10496262;
    % As a symbolic expression, save the expressions 
    % for the complex coefficients below:
    c0 = -0.48209;
    cn = (1i - 1i* e^(-3/2 + 1i*n*pi))/(3*1i + 2*n*pi)+(1i - 1i/e^(1i*n*pi))/(n*pi) + (1i - 1i*e^(-3/2 - 1i*n* pi))/(3*1i - 2*n*pi);
    % Find the complex FS, evaluate the FS up to 
    % second harmonic and save the expression below:
    
    cm2=(1i - 1i* e^(-3/2 + 1i*-2*pi))/(3*1i + 2*-2*pi)+(1i - 1i/e^(1i*-2*pi))/(-2*pi) + (1i - 1i*e^(-3/2 - 1i*-2* pi))/(3*1i - 2*-2*pi);
    em2= e^(2i*pi*-2*1*t);
    cm1=(1i - 1i* e^(-3/2 + 1i*-1*pi))/(3*1i + 2*-1*pi)+(1i - 1i/e^(1i*-1*pi))/(-1*pi) + (1i - 1i*e^(-3/2 - 1i*-1* pi))/(3*1i - 2*-1*pi);
    em1= e^(2i*pi*-1*1*t);
    c1 =(1i - 1i* e^(-3/2 + 1i*1*pi))/(3*1i + 2*1*pi)+(1i - 1i/e^(1i*1*pi))/(1*pi) + (1i - 1i*e^(-3/2 - 1i*1* pi))/(3*1i - 2*1*pi);
    e1= e^(2i*pi*1*1*t);
    c2 =(1i - 1i* e^(-3/2 + 1i*2*pi))/(3*1i + 2*2*pi)+(1i - 1i/e^(1i*2*pi))/(2*pi) + (1i - 1i*e^(-3/2 - 1i*2* pi))/(3*1i - 2*2*pi);
    e2= e^(2i*pi*2*1*t);
    s2t_approx = cm2*em2 + cm1*em1 + c0 + c1*e1 + c2*e2;
end