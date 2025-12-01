%% housekeeping
delete(instrfindall)

%% Stimulator (left coil only)
stimulators = [];
addpath(fullfile('toolboxes','bbci-2e6fe94'), ...
        fullfile('toolboxes','MAGIC-0.1-beta')) % add paths
try
    stimulators.left = magventure('COM3');
    stimulators.left.connect();
catch
    fprintf(2, 'Problem with the MAGIC toolbox, please restart Matlab.\n');
end
%% Initialize BOSS device and spatial filter
bd = bossdevice;
bd.initialize; % not necessary right?

% Stop application if it was alrady running
% if bd.isRunning
%     bd.stop;
% end

bd.start;

%% BD configuration
bd.num_eeg_channels = 64; % Gabor modified this line
% Left-hemisphere Hjorth-like spatial filter (C3 cluster)
spf = zeros(64,1);
spf([5 21 23 25 27]) = [1 -0.25 -0.25 -0.25 -0.25];
bd.spatial_filter_weights = spf;


bd.alpha.offset_samples = 2;           % 2 seems appropriate
bd.alpha.amplitude_min(1) = 0;
bd.alpha.amplitude_max(1) = 1e6;

% ignore other bands (assumes API supports these calls)
bd.theta.ignore;
bd.beta.ignore;

% enforce a minimum inter-trigger interval of 2 seconds
bd.min_inter_trig_interval = 2;

bd.configure_generator_sequence([0 0.001 3 100]); % [time duration port marker]  <-- figure out which port, TS index 100 for BOSS DEVICE

%% Arm stimulator (left)
stimulators.left.arm;
stimulators.left.setAmplitude(15); % can be arbitrary (no actual triggers will be given)

fprintf('Configured stimulator (left coil) and BOSS device.\n');

%% Prepare trigger loop 

n_trials_per_phase = 100; % 100 trials per phase target

n_phases = 3;

% load the random phase values from csv file
random_phases = readtable('data/random_phases.csv');

%% Run trigger loop
fprintf('Entering trigger loop (fixed ITI = %.1f s). Press Ctrl-C to stop.\n', bd.min_inter_trig_interval);

% Phase 1: Trough (pi rad)
fprintf('\nStarting Phase 1: Trough (pi rad)\n');
% initialize the first phase target (pi rad)
bd.alpha.target_phase(1 ) = pi;
bd.alpha.phase_plusminus(1) = pi/40;

% Run Phase dependent stimulation
bd.arm;
while(bd.triggers_remaining > 0)
    pause(0.001); % small pause to prevent overload
end    
bd.disarm;

fprintf('Switching to phase 2: Peak (0 rad).\n');
% configure for second targer
bd.configure_generator_sequence([0 0.001 3 200]); % reconfigure for phase 2 (TS index 200 for NeurOne) 
bd.alpha.target_phase(1) = 0;
bd.triggers_remaining = n_trials_per_phase; % reset trigger count

% Run Phase dependent stimulation
bd.arm;
while(bd.triggers_remaining > 0)
    pause(0.001); % small pause to prevent overload
end
bd.disarm;

fprintf('Switching to phase 3: Random Phase.\n');
bd.configure_generator_sequence([0 0.001 3 300]); % reconfigure for phase 3 (TS index 300 for NeurOne)
random_phase_value = random_phases.phase(total_trials - bd.triggers_remaining + 1);
bd.alpha.target_phase(1) = random_phase_value(1, 1);
bd.triggers_remaining = n_trials_per_phase; % reset trigger count
triggers_before = bd.triggers_remaining;

% Run Phase dependent stimulation
bd.arm;
while(bd.triggers_remaining > 0)
    if triggers_before ~= bd.triggers_remaining
        random_phase_value = random_phases.phase(total_trials - bd.triggers_remaining + 1, 1);
        bd.alpha.target_phase(1) = random_phase_value(1, 1);
        bd.alpha.phase_tolerance = random_phase_value(1, 1)/40;
        triggers_before = bd.triggers_remaining;
    else
        pause(0.001); % small pause to prevent overload
    end
end
bd.disarm;

fprintf('All phases complete. Stopping BOSS device.\n');
bd.stop;