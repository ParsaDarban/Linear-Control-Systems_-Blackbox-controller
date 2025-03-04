clear
clc
close all

%% Unkown system
Unkown_system_output = load('Unkown_system_output.mat', 'Unkown_system_output');

time = Unkown_system_output.Unkown_system_output.time;
values = Unkown_system_output.Unkown_system_output.signals.values;

figure;
plot(time, values, 'b-');
xlabel('Time (s)');
ylabel('Signal Amplitude');
grid on;
hold on

%% First order Estimate
T11 = 4.7;
num1 = 12;
den11 = [T11 1];
G11_s = tf(num1 , den11);
step(G11_s,'g')

T12 = 3.6;
num1 = 12;
den12 = [T12 1];
G12_s = tf(num1 , den12);
step(G12_s,'r')
legend('Unkown system output' , 'First order T=4.7' , 'First order T=3.6')
title('First order Signal');
grid on
hold off

%% Second order estimate
T21 = 2;
T22 = 2.5;
G21_s = tf(12,[T21 1]) * tf(1, [T22 1]);
figure
plot(time, values, 'b-');
hold on 
step(G21_s,'g')

% G22_s = tf(12*1 , [1 , 3,1]);
% step(G22_s,'r')
% 
% G23_s = tf(12*4 , [1 , 12,4]);
% step(G23_s,'k')

G22_s = tf(12*0.3025 , [1 , 1.65,0.3025]);
step(G22_s,'r')

G23_s = tf(12*1.5 , [1 , 7.32,1.5]);
step(G23_s,'k')
grid on

legend('Unkown system output' , 'Second order \zeta=1' , 'Second order \zeta=1.5' , 'Second order \zeta=3')
title('Second order Signal');

%% Estimate base on bode and sine input
w = [0.01, 0.02, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100];
amp = [12, 12, 11.52, 10.39, 5.831, 2.155, 0.5061, 4.78e-2, 7.639e-3, 1.45e-3, 2.055e-4, 5e-5]; 
phase = [0, -5, -25, -48, -100, -150, -180, -210, -230, -240, -270, -270]; 

gain_dB = 20 * log10(amp);
figure;
subplot(2, 1, 1); 
semilogx(w, gain_dB); 
grid on;
xlabel('Frequency (rad/sec)');
ylabel('Magnitude (dB)');
title('Bode Plot - Magnitude');

subplot(2, 1, 2);
semilogx(w, phase);
grid on;
xlabel('Frequency (rad/sec)');
ylabel('Phase (degrees)');
title('Bode Plot - Phase');

% cftool output
ws1 = 1.971;
ws2 = 1.975;  
ws3 = 0.5455;  

G3_s1 = tf(1,[ws1 1]);
G3_s2 = tf(1,[ws2 1]);
G3_s3 = tf(1,[ws3 1]);

G3_s = 12*G3_s1*G3_s2*G3_s3;

% a = 2.097;
% b = 93.71;
% c = 1.576;
% d = 93.49;
% 
% G31_s1 = tf(1,[a 1]);
% G31_s2 = tf(1,[b 1]);
% G31_s3 = tf(1,[c 1]);
% G31_s4 = tf(1,[d 1]);
% G31_s = 12*G31_s1*G31_s2*G31_s3*G31_s4;
% 
% a2 = 1.895e+05;
% b2 = 1.895e+05;  
% 
% G32_s1 = tf(1,[a2 1]);
% G32_s2 = tf(1,[b2 1]);
% G32_s = 12*G32_s1*G32_s2;

figure
plot(time, values, 'b-');
hold on 
step(G3_s,'g')
%step(G32_s,'r')
grid on
legend('Unkown system output' , 'Frequency based') %, 'Frequency based1')
title('Frequency based');

%% Initial G(s)
G1_s = G12_s;
G2_s = G21_s;
G3_s = G3_s;

figure
plot(time, values,'b');
hold on 
step(G1_s,'g')
step(G2_s,'r')
step(G3_s,'y')
grid on
title('All Step Response In one Plot')
legend('Unkown system output','First order','Second order','Frequency Based')

%% Root Locus of Negative feedback
figure;
rlocus(G1_s); 
title('Root Locus with Unity Negative Feedback First Order');
xlabel('Real Axis');
ylabel('Imaginary Axis');

figure
rlocus(G2_s); 
title('Root Locus with Unity Negative Feedback Seconde Order');
xlabel('Real Axis');
ylabel('Imaginary Axis');

figure
rlocus(G3_s); 
title('Root Locus with Unity Negative Feedback Frequency based');
xlabel('Real Axis');
ylabel('Imaginary Axis');


%% Bode
figure;
bode(G1_s); 
grid on
title('Bode First Order');

figure
bode(G2_s); 
title('Bode Seconde Order');
grid on

figure
bode(G3_s); 
title('Bode Frequency based');
grid on

[GM1, PM1, Wcg1, Wcp1] = margin(G1_s);  
disp(['Gain Margin (GM1): ', num2str(GM1), ' (unitless)']);
disp(['Phase Margin (PM1): ', num2str(PM1), ' degrees']);
disp(['Gain Crossover Frequency (Wcg1): ', num2str(Wcg1), ' rad/s']);
disp(['Phase Crossover Frequency (Wcp1): ', num2str(Wcp1), ' rad/s']);

[GM2, PM2, Wcg2, Wcp2] = margin(G2_s);  
disp(['Gain Margin (GM2): ', num2str(GM2), ' (unitless)']);
disp(['Phase Margin (PM2): ', num2str(PM2), ' degrees']);
disp(['Gain Crossover Frequency (Wcg2): ', num2str(Wcg2), ' rad/s']);
disp(['Phase Crossover Frequency (Wcp2): ', num2str(Wcp2), ' rad/s']);

[GM3, PM3, Wcg3, Wcp3] = margin(G3_s);  
disp(['Gain Margin (GM3): ', num2str(GM3), ' (unitless)']);
disp(['Phase Margin (PM3): ', num2str(PM3), ' degrees']);
disp(['Gain Crossover Frequency (Wcg3): ', num2str(Wcg3), ' rad/s']);
disp(['Phase Crossover Frequency (Wcp3): ', num2str(Wcp3), ' rad/s']);

%% Nyquist
figure
nyquist(G1_s);
grid on
title('Nyquist First Order');

figure
nyquist(G2_s);
grid on
title('Nyquist Second Order');

figure
nyquist(G3_s);
grid on
title('Nyquist Frequency based');

%% Controller for G1(s) 
% kc1 = 0.06;
% GC1_s = tf([kc1 kc1*0.87] , [1 0]);
kc1 = 0.1;
GC1_s = tf([kc1 kc1*0.69] , [1 0]);
GD1 = GC1_s*G1_s;
G1c_s = (GD1)/(1+(GD1));
figure
step(G1c_s)
title('Step response with negative feedback for first order')
grid on
info1 = stepinfo(G1c_s);
fprintf('Overshoot = %.3f\n',info1.Overshoot)
fprintf('ts = %.3f\n',info1.SettlingTime)

%% Controller for G2(s) 
% kc2 = 0.33/12;
% GC2_s = tf([kc2 kc2*0.13] , [1 0]);
kc2 = 0.13;
GC2_s = tf([kc2 kc2*0.3] , [1 0]);
% Another choise 
% kc2 = 0.26/1.1;
% GC2_s = tf([kc2 kc2*0.195] , [1 0]);
GD2 = GC2_s*G2_s;
G2c_s = (GD2)/(1+(GD2));
figure
step(G2c_s)
title('Step response with negative feedback')
grid on
info2 = stepinfo(G2c_s);
fprintf('Overshoot = %.3f\n',info2.Overshoot)
fprintf('ts = %.3f\n',info2.SettlingTime)

%% Controller for G3(s)
kc3 = 0.15; 
GC3_s = tf([kc3 kc3*0.2] , [1 0]);
% kc3 = 0.47/12; 
% GC3_s = tf([kc3 kc3*0.015] , [1 0]);
GD3 = GC3_s*G3_s;
G3c_s = (GD3)/(1+(GD3));
figure
step(G3c_s)
title('Step response with negative feedback')
grid on
info3 = stepinfo(G3c_s);
fprintf('Overshoot = %.3f\n',info3.Overshoot)
fprintf('ts = %.3f\n',info3.SettlingTime)

%% All plot in one figure
figure
step(G1c_s)
hold on
step(G2c_s)
step(G3c_s)
grid on
title('All controlled system in one figure')
legend('First order','Second order','Bode estimator')

%% Root locus of controlled system
figure;
rlocus(GD1); 
title('Root Locus with Unity Negative Feedback First Order');
xlabel('Real Axis');
ylabel('Imaginary Axis');

figure
rlocus(GD2); 
title('Root Locus with Unity Negative Feedback Seconde Order');
xlabel('Real Axis');
ylabel('Imaginary Axis');

figure
rlocus(GD3); 
title('Root Locus with Unity Negative Feedback Bode estimator');
xlabel('Real Axis');
ylabel('Imaginary Axis');

%% Bode of controlled system
figure;
bode(GD1); 
grid on
title('Bode First Order');

figure
bode(GD2); 
title('Bode Seconde Order');
grid on

figure
bode(GD3); 
title('Bode for Bode estimator');
grid on

%% Nyquist of controlled system
figure
nyquist(GD1);
grid on
title('Nyquist First Order');

figure
nyquist(GD2);
grid on
title('Nyquist Second Order');

figure
nyquist(GD3);
grid on
title('Nyquist Bode estimator');

%% Controller of the BlackBox
FirstOrderController = load('FirstOrderController.mat', 'FirstOrderController');
timeFirst = FirstOrderController.FirstOrderController.time;
valuefirst = FirstOrderController.FirstOrderController.signals.values;

SecondOrderController = load('SecondOrderController.mat', 'SecondOrderController');
timeSecond = SecondOrderController.SecondOrderController.time;
valueSecond = SecondOrderController.SecondOrderController.signals.values;

BodeBaseController = load('BodeBaseController.mat', 'BodeBaseController');
timeBode = BodeBaseController.BodeBaseController.time;
valueBode = BodeBaseController.BodeBaseController.signals.values;

figure;
plot(timeFirst, valuefirst);
hold on
plot(timeSecond, valueSecond);
plot(timeBode, valueBode);
xlabel('Time (s)');
ylabel('Signal Amplitude');
grid on
legend('FirstOrderController','SecondOrderController','BodeBaseController')

%% Frequency Controller
zeta = 0.591;
PM = 100 * zeta;
figure
margin(G2_s)
grid on
k = 0.125;
G22_s = G2_s;
w = 0.31;
Gfc_s = tf([k k*w],[1 0]);
GDf = Gfc_s*G22_s;
Gfc = (GDf)/(1+(GDf));
figure
step(Gfc)
grid on
stepinfo(Gfc)

figure
bode(G2_s);
legend('system bode', 'feedback bode')
hold on
bode(GDf)
grid on
legend('system bode', 'feedback bode')

[GM2, PM2, Wcg2, Wcp2] = margin(G2_s);  
disp(['Gain Margin (GM2): ', num2str(GM2), ' (unitless)']);
disp(['Phase Margin (PM2): ', num2str(PM2), ' degrees']);
disp(['Gain Crossover Frequency (Wcg2): ', num2str(Wcg2), ' rad/s']);
disp(['Phase Crossover Frequency (Wcp2): ', num2str(Wcp2), ' rad/s']);

%% Black Box frequency Controlled
FrequencyBasedController = load('FrequencyBasedController.mat', 'FrequencyBasedController');
timeFrequency = FrequencyBasedController.FrequencyBasedController.time;
valueFrequency = FrequencyBasedController.FrequencyBasedController.signals.values;
figure;
plot(timeFrequency, valueFrequency);
grid on
title('Frequnecy Based Controlled system')
xlabel('Time (s)');
ylabel('Signal Amplitude');

%% Function for this project that may use
% function k = stable(G)
%     K_values = linspace(0, 1000, 10000000);
%     
%     poles = rlocus(G, K_values);
%     
%     stable_K = K_values(all(real(poles) < 0, 1));
%     
%     if isempty(stable_K)
%         disp('No stable range for K.');
%         k = []; 
%     else
%         k = [min(stable_K), max(stable_K)]; 
%         disp(['Stable K range: ', num2str(k(1)), ' to ', num2str(k(2))]);
%     end
% end
