%% housekeeping
delete(instrfindall)
clc;clear;

% Load settings from JSON file
settings = jsondecode(fileread('C:\Users\Eric James McDermott\Desktop\LAVA_GitHub\code\settings.json'));

%% Initialize BOSS device and spatial filter
bd = bossdevice;
pause(1);
% Stop application if it was alrady running
if bd.isRunning
    bd.stop;
end
pause(1);
bd.start;

%% BD configuration
bd.num_eeg_channels = 64; 

% Left-hemisphere Hjorth-like spatial filter (C3 cluster)
spf = zeros(64,1);
spf([5 21 23 25 27]) = [1 -0.25 -0.25 -0.25 -0.25];
bd.spatial_filter_weights = spf;
coeffs = load([settings.data_save '\bpfilter_' settings.sub_id '.mat'], 'coefficients').coefficients;
bd.alpha.bpf_fir_coeffs = coeffs;

bd.alpha.offset_samples = 2;    % 2 seems appropriate
bd.alpha.amplitude_min(1) = 0;
bd.alpha.amplitude_max(1) = 1e6;

% ignore other bands (assumes API supports these calls)
bd.theta.ignore;
bd.beta.ignore;

% enforce a minimum inter-trigger interval of 2 seconds
bd.min_inter_trig_interval = 2;

%% Prepare trigger loop 

% load the random phase values from csv file
phase_targets = readtable([settings.data_save '\phase_schedule.csv']);

% get length of random phase array
n_trials = height(phase_targets);

%% Run trigger loop
fprintf('Entering trigger loop (fixed ITI = %.1f s). Press Ctrl-C to stop.\n', bd.min_inter_trig_interval);

bd.alpha.phase_plusminus(1) = pi/40;

% safeguard 
bd.triggers_remaining = 0;

% Arm the BOSS Device
bd.arm;

% waiting for Neurosimo to sent a udp package to trigger BOSS loop
matlab_udp_listener(5555);

% Configure phase trials
bd.triggers_remaining = n_trials; % reset trigger count
triggers_before = bd.triggers_remaining;

random_phase_value = phase_targets.Var1(n_trials - bd.triggers_remaining + 1, 1);
bd.alpha.phase_target(1) = random_phase_value(1, 1);

% Run Phase dependent stimulation
bd.arm;
while(bd.triggers_remaining > 0)
    if triggers_before ~= bd.triggers_remaining && bd.triggers_remaining > 0 % additional safeguard to avoid index error in table

        phase_value = phase_targets.Var1(n_trials - bd.triggers_remaining + 1, 1);
        
        % Map phase target to marker value
        if phase_value == 0
            marker_val = 10;
        elseif abs(phase_value - pi) < 1e-4
            marker_val = 20;
        else
            marker_val = 30;
        end
        
        bd.configure_generator_sequence([0 0.001 marker_val]);
        bd.alpha.phase_target(1) = phase_value;
        triggers_before = bd.triggers_remaining;
    else
        pause(0.001); % small pause to prevent overload
    end
end

bd.disarm;
fprintf('All trials complete. Stopping BOSS device.\n');

%% End experiment
bd.stop;