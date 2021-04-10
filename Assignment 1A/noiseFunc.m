function [t, s2_hinf, sid] = noiseFunc(sid);

    
    % Enter you student ID below:
    sid = 10496262;
    % Save the appropriate outputs below as defined above:
    %samples = 100;
    t = linspace(-1/2 , 1/2, 100 + 1); t(end) = [];
    s2t = zeros(1,100);
    s2t(t >= -1/2 & t < 0) = exp(3*t(t >= -1/2 & t < 0));
    s2t(t >= 0 & t < 1/2) = exp(-3*t(t >= 0 & t < 1/2)) - 2; % I didnt suppress the output on my other attempts
    s2_hinf = repmat(s2t,1,5);
    t = linspace(-1/2 , 4.5, length(s2t)+1); t(end) = [];

    
    %Plot function
    figure();
    plot(t, s2_hinf, 'b', 'LineWidth', 2);

    
end