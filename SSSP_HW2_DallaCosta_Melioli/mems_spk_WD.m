%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SSSP HOMEWORK - SASP MODULE 2, A.Y 24/25 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DALLA COSTA MATTIA - MELIOLI ANNA CHIARA  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
Z2 = R1;
Z5 = Ts/(2*C1);
Z8= R2;
Z9 = 2*L1/Ts;
Z10 = Ts/(2*C2);
Z13 = Ts/(2*C3);
Z16 = 2*L2/Ts;
Z19 = Ts/(2*C4);
Z22 = 2*L3/Ts;
Z23 = R3;

Z21 = Z22 + Z23;
Z20 = Z21;
Z18 = Z19*Z20/(Z19+Z20);
Z17 = Z18;
Z15 = Z17 + Z16;
Z14 = Z15;
Z12 = Z14*Z13/(Z14 + Z13);
Z11 = Z12;
Zint1 = Z11 + Z10;
Zint2 = Zint1;
Zint3 = Zint2 + Z9;
Zint4 = Zint3;
Z7 = Z8 + Zint4;
%Z7 = Z8 + Z9 + Z10;
Z6 = Z7;
Z4 = Z5*Z6/(Z5 + Z6);
Z3 = Z4;
Z1 = Z2 + Z3;


%% Computing Scattering Matrices
B1 = [1, 1, 1];    
% B3 = [1, 1, 1, 1, 1];
B31 = [1, 1, 1];
B32 = [1, 1, 1];
B33 = [1, 1, 1];
B5 = [-1, 1, 1];
B7 = [-1, 1, 1];
Q = [1, 1, 1];
tic
Zser1 = diag([Z1, Z2, Z3]);
S1 = eye(3) - 2*Zser1 * B1' * inv(B1 * Zser1 * B1') * B1;

Zpar2 = diag([Z4, Z5, Z6]);
S2 = 2 * Q' * inv(Q * inv(Zpar2) * Q') * Q * inv(Zpar2) - eye(3);

% Zser3 = diag([Z7, Z8, Z9, Z10, Z11]);
% S3 = eye(5) - 2*Zser3 * B3' * inv(B3 * Zser3 * B3') * B3;
Zser31 = diag([Z7,Z8,Zint4]);
S31 = eye(3) - 2*Zser31 * B31' * inv(B31 * Zser31 * B31') * B31;
Zser32 = diag([Zint3,Z9,Zint2]);
S32 = eye(3) - 2*Zser32 * B32' * inv(B32 * Zser32 * B32') * B32;
Zser33 = diag([Zint1,Z10,Z11]);
S33 = eye(3) - 2*Zser33 * B33' * inv(B33 * Zser33 * B33') * B33;

Zpar4 = diag([Z12, Z13, Z14]);
S4 = 2 * Q' * inv(Q * inv(Zpar4) * Q') * Q * inv(Zpar4) - eye(3);

Zser5 = diag([Z15, Z16, Z17]);
S5 = eye(3) - 2*Zser5 * B5' * inv(B5 * Zser5 * B5') * B5;

Zpar6 = diag([Z18, Z19, Z20]);
S6 = 2 * Q' * inv(Q * inv(Zpar6) * Q') * Q * inv(Zpar6) - eye(3);

Zser7 = diag([Z21, Z22, Z23]);
S7 = eye(3) - 2*Zser7 * B7' * inv(B7 * Zser7 * B7') * B7;

%% Initialization of Waves
a = zeros(23,1);
b = zeros(23,1);
aint = zeros(4,1);
bint = zeros(4,1);

%% Initialization of Output Signals
Fout = zeros(1, length(t));

%% Simulation Algorithm
for n = 1:length(Fin)
    % Relations
    a(22) = -b(22);    % L3
    a(19) = b(19);     % C4
    a(16) = -b(16);    % L2
    a(13) = b(13);     % C3
    a(10) = b(10);     % C2
    a(9) = -b(9);      % L1
    a(5) = b(5);       % C1

    % Forward Scan
    b(21) = S7(1,:) * [0; a(22); a(23)];
    a(20) = b(21);
    b(18) = S6(1,:) * [0; a(19); a(20)];
    a(17) = b(18);
    b(15) = S5(1,:) * [0; a(16); a(17)];
    a(14) = b(15);
    b(12) = S4(1,:) * [0; a(13); a(14)];
    a(11) = b(12);
    bint(1) = S33(1,:) * [0; a(10); a(11)];
    aint(2) = bint(1);
    bint(3) = S32(1,:) * [0; a(9); aint(2)];
    aint(4) = bint(3);
    b(7) = S31(1,:) * [0; a(8); aint(4)];
    % b(7) = S3(1,:) * [0; a(8); a(9); a(10); a(11)];

    a(6) = b(7);
    b(4) = S2(1,:) * [0; a(5); a(6)];
    a(3) = b(4);
    b(1) = S1(1,:) * [0; a(2); a(3)];
    
    % Local Root Scattering
    aroot = b(1);
    broot = 2*Fin(n) - aroot;
    a(1) = broot;

    % Backward Scan
    b(2) = S1(2,:) * [a(1); a(2); a(3)];
    b(3) = S1(3,:) * [a(1); a(2); a(3)];
    a(4) = b(3);

    b(5) = S2(2,:) * [a(4); a(5); a(6)];
    b(6) = S2(3,:) * [a(4); a(5); a(6)];
    a(7) = b(6);
    
    b(8) = S31(2,:) * [a(7); a(8); aint(4)];
    bint(4) = S31(3,:) * [a(7); a(8); aint(4)];
    aint(3) = bint(4);
    b(9) = S32(2,:) * [aint(3); a(9); aint(2)];
    bint(2) = S32(3,:) * [aint(3); a(9); aint(2)];
    aint(1) = bint(2);
    b(10) = S33(2,:) * [aint(1); a(10); a(11)];
    b(11) = S33(3,:) * [aint(1); a(10); a(11)];

    %b(8) = S3(2,:) * [a(7); a(8); a(9); a(10); a(11)];
    %b(9) = S3(3,:) * [a(7); a(8); a(9); a(10); a(11)];
    %b(10) = S3(4,:) * [a(7); a(8); a(9); a(10); a(11)];
    %b(11) = S3(5,:) * [a(7); a(8); a(9); a(10); a(11)];
    a(12) = b(11);

    b(13) = S4(2,:) * [a(12); a(13); a(14)];
    b(14) = S4(3,:) * [a(12); a(13); a(14)];
    a(15) = b(14);

    b(16) = S5(2,:) * [a(15); a(16); a(17)];
    b(17) = S5(3,:) * [a(15); a(16); a(17)];
    a(18) = b(17);

    b(19) = S6(2,:) * [a(18); a(19); a(20)];
    b(20) = S6(3,:) * [a(18); a(19); a(20)];
    a(21) = b(20);

    b(22) = S7(2,:) * [a(21); a(22); a(23)];
    b(23) = S7(3,:) * [a(21); a(22); a(23)];

    % Read Output
    Fout(n) = (b(23) + a(23)) / 2;
end

%% Output Plots

% Computing acoustic pressure
pout = Fout ./ Seff;
toc
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
