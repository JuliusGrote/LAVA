%% housekeeping
delete(instrfindall)
clc;clear;
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
coeffs = load('data/filter_coeffs.mat', 'coeffs').coeffs;
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

n_trials_per_phase = 200; % 100 trials per phase target

n_phases = 3;

% load the random phase values from csv file
random_phases = readtable('data/random_phases.csv');

%% Run trigger loop
fprintf('Entering trigger loop (fixed ITI = %.1f s). Press Ctrl-C to stop.\n', bd.min_inter_trig_interval);

% Phase 1: Trough (pi rad)
fprintf('\nStarting Phase 1: Trough (pi rad)\n');
% initialize the first phase target (pi rad)
bd.configure_generator_sequence([0 0.001 3 10]); % [time duration port marker], TS index 100 for BOSS DEVICE phase target "pi"
bd.alpha.phase_target(1) = pi;
bd.alpha.phase_plusminus(1) = pi/40;

% safeguard 
bd.triggers_remaining = 0;

% Arm the BOSS Device
bd.arm;

% waiting for Neurosimo to sent a udp package to trigger BOSS loop
matlab_udp_listener(5555);
bd.triggers_remaining = n_trials_per_phase;

while(bd.triggers_remaining > 0)
    pause(0.001); % small pause to prevent overload
end    
bd.disarm;

fprintf('Switching to phase 2: Peak (0 rad).\n');
% configure for second targer
bd.configure_generator_sequence([0 0.001 3 20]); % reconfigure for phase 2 (TS index 200 for phase target "0") 
bd.alpha.phase_target(1) = 0;
bd.triggers_remaining = n_trials_per_phase; % reset trigger count

% Run Phase dependent stimulation
bd.arm;
while(bd.triggers_remaining > 0)
    pause(0.001); % small pause to prevent overload
end
bd.disarm;

% Configure random phase trials
bd.triggers_remaining = n_trials_per_phase; % reset trigger count
triggers_before = bd.triggers_remaining;
fprintf('Switching to phase 3: Random Phase.\n');
bd.configure_generator_sequence([0 0.001 3 30]); % reconfigure for phase 3 (TS index 300 for random phase target)
random_phase_value = random_phases.Var1(n_trials_per_phase - bd.triggers_remaining + 1, 1);
bd.alpha.phase_target(1) = random_phase_value(1, 1);


% Run Phase dependent stimulation
bd.arm;
while(bd.triggers_remaining > 0)
    if triggers_before ~= bd.triggers_remaining && bd.triggers_remaining > 0 % additional safeguard to avoid index error in table

        random_phase_value = random_phases.Var1(n_trials_per_phase - bd.triggers_remaining + 1, 1);
        bd.alpha.phase_target(1) = random_phase_value;
        triggers_before = bd.triggers_remaining;
    else
        pause(0.001); % small pause to prevent overload
    end
end
bd.disarm;

fprintf('All phases complete. Stopping BOSS device.\n');

%% End experiment
bd.stop;