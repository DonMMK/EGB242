t0_3 = 0;
tend_3 = 20;
T3 = 5;
x = linspace(t0_3 , t0_3 + T3, 100 + 1); x(end) = [];
t3 = repmat(x, [1 4]);
func = @(t3) (A*t3)/4.*(0 >= t3 & t3 < 2.5) + (-A*t3 + 5*A)/4 .*(2.5 >= t3 & t3 < 5);
%s3_t = func1(t3) + func2(t3);
t3 = linspace(t0_3 , t0_3 + 4*T3, 400 + 1); t3(end) = [];
plot(t3,func(t3), 'r', 'LineWidth', 2)


% set(0,'DefaultFigureWindowStyle','docked');
% Coeficients
% A = 5; B = 3;
% time = 20;
% T = time/4;
% samples = 44100; %per second
% fivesec = samples/5;%for 5 seconds
% 
% x = linspace(t0 , t0 + T, 100 + 1); x(end) = [];
% t = repmat(x, [1 5]);
%  
% sampleperiod = samples * T;
% additive = noiseSound(1:sampleperiod);
%  
% t = linspace(0, T, sampleperiod + 1);t(end) =[];
%  
% s1_t = t/A; 
% s2_t = exp(-(t-B)/4);
% func1 = @(t) (A*t)/4.*(0 >= t) & (t < 2.5) + 0.0;
% func2 = @(t) exp((-A*t + 5*A)/4) .*(2.5 >= t) & (t < 5) + 0.0;
% s3_t = func1(t) + func2(t);
% 
% s3_t( t >= 0 & t < 2.5 ) = (A*t)/4;
% s3_t( t >= 2.5 & t < 5 ) = (-A*t + 5*A)/4;
%      
% figure()
% plot(t, s1_t, 'r', 'LineWidth', 2);
% hold on
% plot(t, s2_t, 'b', 'LineWidth', 2);
% plot(t, s3_t, 'k', 'LineWidth', 2);
% plot(t, additive, 'm', 'LineWidth', 2);
% grid on
% title('Plot', 'FontSize', 30);
