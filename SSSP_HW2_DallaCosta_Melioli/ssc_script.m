clear; close all; clc

%% Simscape Controlling Script
% Open the Simscape model 'mems_spk_SSC.slx' and run the following commands
% to define the Simscape Input Signal as a MATLAB workspace variable.

% Sampling Frequency
fs = 192e3;

% Sampling Period
Ts = 1/fs;

% Simulation Duration
stop_time = 2;  % [seconds]

% Input Signal

% Time Axis
t = 0:Ts:stop_time;
t = t(1:end-1);

% Signal Amplitude
A = 10;
vin = A * sweeptone(1.99, 0.01, fs, 'SweepFrequencyRange', [500 20000]);

% Define Simscape Input Signal (ssc_input)
ssc_input = [t', vin];