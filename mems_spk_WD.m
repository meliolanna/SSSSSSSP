clear; close all; clc

%% Simscape File for Ground Truth

load('ssc_output.mat')

%% Sampling Frequency
fs = 192e3;

%% Sampling Period
Ts = 1/fs;

%% Simulation Duration
stop_time = 2;  % [seconds]

%% Input Signal

% Time Axis
t = 0:Ts:stop_time;
t = t(1:end-1);

% Signal Amplitude
A = 10;
% vin = A * sweeptone(1.99, 0.01, fs, 'SweepFrequencyRange', [500 20000]);
load('vin.mat');

%% Circuit Parameters
% Piezoelectric MEMS loudspeaker in Free-field conditions

% Transformer Ratio
alpha = 3.7e-4;
Seff = 2.0e-5;

% Electrical Domain
Re = 4;
Cp = 2.4e-8;

% Mechanical Domain
Rm = 9.7e-3;
Mm = 1e-6;
Cm = 2.2e-3;

% Acoustic Domain
Cbc = 3.6e-13; 
Ltube1 = 1e2;
Ltube2 = 1e2;
Ctube = 6.5e-13;
Rac = 5e6;

%% Removing Ideal Transformers (Mechanical Domain only)

gamma1 = alpha^-1;
gamma2 = Seff;

% Resistive Elements
R1 = Re/gamma1^2;
R2 = Rm;
R3 = Rac*gamma2^2;


% Dynamic Elements
% Capacitors
C1 = Cp*gamma1^2;
C2 = Cm;
C3 = Cbc/gamma2^2;
C4 = Ctube/gamma2^2;


% Inductors
L1 = Mm;
L2 = Ltube1*gamma2^2;
L3 = Ltube2*gamma2^2;

% Source
Fin = alpha * vin;

%% Setting of Free Parameters (Adaptation Conditions)
Z1 = R3;
Z2 = 2*L3/Ts;
Z5 = Ts/(2*C4);
Z8 = 2*L2/Ts;
Z11 = Ts/(2*C3);
Z14 = Ts/(2*C2);
Z15 = R2;
Z16 = 2*L1/Ts;
Z19 = Ts/(2*C1);
Z22 = R1;

Z3 = Z1 + Z2;
Z4 = Z3;
Z6 = Z4*Z5 / (Z4+Z5);
Z7 = Z6;
Z9 = Z7 + Z8;
Z10 = Z9;
Z12 = Z10*Z11 / (Z10 + Z11);
Z13 = Z12;
Z17 = Z14 + Z15 + Z16;
Z18 = Z17;
Z20 = Z18*Z17/(Z18+Z17);
Z21 = Z20;
Z23 = Z21+Z22;

%% Computing Scattering Matrices
B1 = [-1, 1, 1];
B3 = B1;
B5 = [-1, 1, 1, 1, 1];
B7 = [1, 1, 1];
Q2 = [1, 1, 1];
Q4 = Q2;
Q6 = Q2;

Zser1 = diag([Z1, Z2, Z3]);
S1 = eye(3) - 2*Zser1*(B1.')*inv(B1*Zser1*(B1.'))*B1;

Zpar2 = diag([Z4, Z5, Z6]);
S2 = 2*(Q2.')*inv(Q2*inv(Zpar2)*(Q2.'))*Q2*inv(Zpar2) - eye(3);

Zser3 = diag([Z7, Z8, Z9]);
S3 = eye(3) - 2*Zser3*B3'*inv(B3*Zser3*B3')*B3;

Zpar4 = diag([Z10, Z11, Z12]);
S4 = 2*(Q4.')*inv(Q4*inv(Zpar4)*(Q4.'))*Q4*inv(Zpar4) - eye(3);

Zser5 = diag([Z13, Z14, Z15, Z16, Z17]);
S5 = eye(5) - 2*Zser5*B5'*inv(B5*Zser5*B5')*B5;

Zpar6 = diag([Z18, Z19, Z20]);
S6 = 2*(Q6.')*inv(Q6*inv(Zpar6)*(Q6.'))*Q6*inv(Zpar6) - eye(3);

Zser7 = diag([Z21, Z22, Z23]);
S7 = eye(3) - 2*Zser7*B7'*inv(B7*Zser7*B7')*B7;

%% Initialization of Waves
a = zeros(23,1);
b = zeros(23,1);

%% Initialization of Output Signals
Fout = zeros(1, length(t));

%% Simulation Algorithm

for n = 1 : length(Fin)
    % Relations
    b(2) = - a(2);
    b(5) = a(5);
    b(8) = - a(8);
    b(11) = a(11);
    b(14) = a(14);
    b(16) = a(16);
    b(19) = a(19);

    % Forward Scan
    b(3) = Zser1(3,:)*[a(1),a(2),0].';
    a(4) = b(3);    
    b(6) = Zpar2(3,:)*[a(4),a(5),0].';
    a(7) = b(6);   
    b(7) = Zser3(3,:)*[a(7),a(8),0].';
    a(8) = b(7);   
    b(10) = Zpar4(3,:)*[a(10),a(11),0].';
    a(11) = b(10);   
    b(17) = Zser5(5,:)*[a(13),a(14),a(15),a(16),0].';
    a(18) = b(17);
    b(20) = Zpar6(3,:)*[a(18),a(19),0].';
    a(21) = b(20);
    b(23) = Zser7(3,:)*[a(21),a(22),0].';

    % Local Root Scattering
    aroot = b(23);
    broot = 2*Fin(n) - aroot;
    a(23) = broot;

    % Backward Scan
    b(21) = Zser7(1,:)*a(21:23);
    b(22) = Zser7(2,:)*a(21:23);
    a(20) = b(21);
    b(19) = Zpar6(2,:)*a(18:20);
    b(18) = Zpar6(1,:)*a(18:20);
    a(17) = b(18);
    b(16) = Zser5(4,:)*a(13:17);
    b(15) = Zser5(3,:)*a(13:17);
    b(14) = Zser5(2,:)*a(13:17);
    b(13) = Zser5(1,:)*a(13:17);
    a(12) = b(13);
    b(11) = Zpar4(2,:)*a(10:12);
    b(10) = Zpar4(1,:)*a(10:12);
    a(9) = b(10);
    b(8) = Zser3(2,:)*a(7:9);
    b(7) = Zser3(1,:)*a(7:9);
    a(6) = b(7);
    b(5) = Zpar2(2,:)*a(4:6);
    b(4) = Zpar2(1,:)*a(4:6);
    a(3) = b(4);
    b(2) = Zser1(2,:)*a(1:3);
    b(1) = Zser1(1,:)*a(1:3);

    % Read Output
    Fout(n) = (b(1)+a(1))/2;
end

%% Output Plots

% Computing acoustic pressure
pout = Fout ./ Seff;

% Time Domain Plots
figure
set(gcf, 'Color', 'w');
plot(t, pout, 'b', 'Linewidth', 2);
hold on
plot(gt(1, :), gt(2, :), 'r--', 'Linewidth', 2);
grid on;
xlabel('Time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$p_{\mathrm{out}}$ [Pa]','Fontsize',16,'interpreter','latex');
legend('WDF','SSC','Fontsize',16,'interpreter','latex');
title('Output Pressure - Time Domain','Fontsize',18,'interpreter','latex');

% Frequency domain Plots
nfft = 2^20;
res = fs/nfft;
f = (0:nfft/2-1) * res;

ir = impzest(vin/A, pout');
ir_gt = impzest(vin/A, gt(2, 1:end-1)');

tf = fft(ir, nfft);
tf_gt = fft(ir_gt, nfft);

abs_tf = abs(tf(1:nfft/2));
abs_tf_gt = abs(tf_gt(1:nfft/2));

figure
set(gcf, 'Color', 'w');
semilogx(f, 20*log10(abs_tf/2e-5), 'b', 'Linewidth', 2);
hold on
semilogx(f, 20*log10(abs_tf_gt/2e-5), 'r--', 'Linewidth', 2);
grid on;
xlim([500, 20000])
xlabel('Frequency [Hz]','Fontsize',16,'interpreter','latex');
ylabel('$\mathrm{SPL}\,[\mathrm{dB}_\mathrm{SPL}]$','Fontsize',16,'interpreter','latex');
legend('WDF','SSC','Fontsize',16,'interpreter','latex');
title('Output Sound Pressure Level - Frequency Domain','Fontsize',16,'interpreter','latex');

%% Error Plots

figure
set(gcf, 'Color', 'w');
hold on;
plot(t, pout - gt(2, 1:end-1), 'k', 'Linewidth', 2);
grid on;
xlabel('Time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$\mathcal{E}_{\mathrm{out}}$ [Pa]','Fontsize',16,'interpreter','latex');
title('Error Signal','Fontsize',16,'interpreter','latex');

%% Compute Mean Squared Error (MSE)

mse = mean((pout - gt(2, 1:end-1)).^2);
disp('MSE = ')
disp(mse)
