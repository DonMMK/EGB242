function [c0, cn, s2t_approx, sid] = compFS(sid)
    
    % Enter your student ID as sid:
    sid = 10496262;
    % As a symbolic expression, save the expressions 
    % for the complex coefficients below:
    c0 = -0.48209;
    %cn = -(1i *(-1 + exp^(-1i*pi* n)))/(pi* n) - (exp^(3/2) - exp^(1i*pi* n))/(exp^(3/2)*(-3 + 2*1i*pi* n)) + (1 - exp^(-3/2 - 1i*pi* n))/(3 + 2 *1i*pi*n);
    cn = (1i - 1i* exp^(-3/2 + 1i*n*pi))/(3*1i + 2*n*pi)+(1i - 1i/exp(1i*n*pi))/(n*pi) + (1i - 1i*exp(-3/2 - 1i*n* pi))/(3*1i - 2*n*pi);
    % Find the complex FS, evaluate the FS up to 
    % second harmonic and save the expression below:
    s2t_approx = (-(1i *(-1 + exp^(-1i*pi* -2)))/(pi* -2) - (exp^(3/2) - exp^(1i*pi* -2))/(exp^(3/2)*(-3 + 2*1i*pi* -2)) + (1 - exp^(-3/2 - 1i*pi* -2))/(3 + 2 *1i*pi*-2))*exp(1i*2*pi*-2*t)  +          (-(1i *(-1 + exp^(-1i*pi* -1)))/(pi* -1) - (exp^(3/2) - exp^(1i*pi* -1))/(exp^(3/2)*(-3 + 2*1i*pi* -1)) + (1 - exp^(-3/2 - 1i*pi* -1))/(3 + 2 *1i*pi*-1))*exp(1i*2*pi*-1*t)    + (-(1i *(-1 + exp^(-1i*pi* 0)))/(pi* 0) - (exp^(3/2) - exp^(1i*pi* 0))/(exp^(3/2)*(-3 + 2*1i*pi* 0)) + (1 - exp^(-3/2 - 1i*pi* 0))/(3 + 2 *1i*pi*0))*exp(1i*2*pi*0*t)  + (-(1i *(-1 + exp^(-1i*pi* 1)))/(pi* 1) - (exp^(3/2) - exp^(1i*pi* 1))/(exp^(3/2)*(-3 + 2*1i*pi* 1)) + (1 - exp^(-3/2 - 1i*pi* 1))/(3 + 2 *1i*pi*1))*exp(1i*2*pi*1*t)  + (-(1i *(-1 + exp^(-1i*pi* 2)))/(pi* 2) - (exp^(3/2) - exp^(1i*pi* 2))/(exp^(3/2)*(-3 + 2*1i*pi* 2)) + (1 - exp^(-3/2 - 1i*pi* 2))/(3 + 2 *1i*pi*2))*exp(1i*2*pi*2*t);  
end
% 
% function [c0, cn, s2t_approx, sid] = compFS(sid)
%     syms n t; e = sym(exp(sym(1))); pi = sym('pi');
%     % Enter your student ID as sid:
%     sid = 10406689;
%     % As a symbolic expression, save the expressions 
%     % for the complex coefficients below:
%     c0 = (1/6)*(-2*e^(-3)-4);
%     cn = (1/(6-2i*pi*n))*(1-e^(-3+pi*n*1i))+((e^(-3-pi*n*1i)-1)/(-6-2i*n*pi)+(e^(-pi*n*1i)-1)/(pi*n*1i));
%     % Find the complex FS, evaluate the FS up to 
%     % second harmonic and save the expression below:
%     a1 = -2;
%     a2 = -1;
%     a3 = 0;
%     a4 = 1;
%     a5 = 2;
%     
%     
%     cnm2 = (1/(6-2i*pi*a1))*(1-e^(-3+pi*a1*1i))+((e^(-3-pi*a1*1i)-1)/(-6-2i*a1*pi)+(e^(-pi*a1*1i)-1)/(pi*a1*1i));
%     cnm1 = (1/(6-2i*pi*a2))*(1-e^(-3+pi*a2*1i))+((e^(-3-pi*a2*1i)-1)/(-6-2i*a2*pi)+(e^(-pi*a2*1i)-1)/(pi*a2*1i));
%  %no cn0 
%     cn1  = (1/(6-2i*pi*a4))*(1-e^(-3+pi*a4*1i))+((e^(-3-pi*a4*1i)-1)/(-6-2i*a4*pi)+(e^(-pi*a4*1i)-1)/(pi*a4*1i));
%     cn2  = (1/(6-2i*pi*a5))*(1-e^(-3+pi*a5*1i))+((e^(-3-pi*a5*1i)-1)/(-6-2i*a5*pi)+(e^(-pi*a5*1i)-1)/(pi*a5*1i));
%     
%     aprox1 = cnm2 * e^(2i*pi*a1*t*(1/3))
%     aprox2 = cnm1 * e^(2i*pi*a2*t*(1/3))
%     %aprox3 = cn0 
%     aprox4 = cn1 * e^(2i*pi*a4*t*(1/3))
%     aprox5 = cn2 * e^(2i*pi*a5*t*(1/3))
%     
%     s2t_approx =  aprox1 +  aprox2 + c0 +  aprox4 + aprox5
%   % st2a = - e^(-3)/3 + e^(-(pi*t*4i)/3)*((e^(-3) - 1)/(- 6 + pi*4i) - (e^(-3) - 1)/(6 + pi*4i)) + e^((pi*t*4i)/3)*((e^(-3) - 1)/(- 6 + pi*4i) - (e^(-3) - 1)/(6 + pi*4i)) - e^(-(pi*t*2i)/3)*((e^(-3) + 1)/(- 6 + pi*2i) - (e^(-3) + 1)/(6 + pi*2i) + 2i/pi) + e^((pi*t*2i)/3)*(- (e^(-3) + 1)/(- 6 + pi*2i) + (e^(-3) + 1)/(6 + pi*2i) + 2i/pi) - 2/3
%  
%   %s2ta = 1/3 + e^(-(pi*t*4i)/3)*((e^(-3) - 1)/(- 6 + pi*4i) - (e^(-3) - 1)/(6 + pi*4i)) + e^((pi*t*4i)/3)*((e^(-3) - 1)/(- 6 + pi*4i) - (e^(-3) - 1)/(6 + pi*4i)) - e^(-(pi*t*2i)/3)*((e^(-3) + 1)/(- 6 + pi*2i) - (e^(-3) + 1)/(6 + pi*2i) + 2i/pi) + e^((pi*t*2i)/3)*(- (e^(-3) + 1)/(- 6 + pi*2i) + (e^(-3) + 1)/(6 + pi*2i) + 2i/pi) - e^(-3)/3
% end