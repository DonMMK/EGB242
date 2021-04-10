function [t, s2_hinf, sid] = noiseFunc(sid);
    % Enter you student ID below:
    sid = 10406689 ;
    % Save the appropriate outputs below as defined above:
    t = linspace(-3/2,3/2,100+1); %t values 
    t(end) = [];
    
    s2 = zeros(1,100);
    s2(t>=-3/2 & t<0) = exp(2*t(t>=-3/2 & t<0));
    s2(t>=0 & t<3/2) = exp(-2*t(t>=0 & t<3/2))-2;
    
    figure
    plot(t,s2);
    
    s2_hinf = repmat(s2,1,5);
    t = linspace(-3/2,27/2,length(s2_hinf)+1); %t values 
    t(end) = [];
    
    figure
    plot(t,s2_hinf);
    

end

%  %% alt
%     % Enter you student ID below:
%     sid = 10496262;
%     % Save the appropriate outputs below as defined above:
%     %samples = 100;
%     t0 = -1/2;
%     tend = 1/2;
%     T = 1;
%     t = linspace(-1/2 , 1/2, 100 + 1); t(end) = [];
%     
%     
%     s2_func = exp(3*t).*(t >= -1/2 & t < 0);
%     s2_func(t >= 0 & t < 1/2) = exp(-3*t(t >= 0 & t < 1/2) - 2); % I didnt suppress the output on my other attempts
%     t = linspace(-1/2 , 5/2, 500 + 1); t(end) = [];
%     s2_hinf = repmat(s2_func, 1, 5);
%     
%     %Plot function
%     figure();
%     plot(t, s2_hinf, 'y', 'LineWidth', 5);
%     hold on
%     plot(t, s2_hinf, 'r:', 'LineWidth', 2);
% 
%     
%  end
% 
% 
% 
