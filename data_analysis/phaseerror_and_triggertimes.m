%% Read NeurOne raw EEG session
data_load = '/home/lisa/projects/LAVA/data/2025-11-21T111154';
data_save = '/home/lisa/projects/LAVA/data';
data_csv = '/home/lisa/projects/LAVA/data/random_phases.csv';
fig_dir = '/home/lisa/projects/LAVA/figures/';

% Define colorblind-friendly palette (high contrast, distinct for all types of colorblindness)
color_NS = [216/255, 90/255, 27/255];         % Orange rgba(216, 90, 27, 1) - distinguishable by all
color_BOSS = [0/255, 114/255, 178/255];       % Blue rgb(0, 114, 178) - safe for colorblind
color_line = [34/255, 34/255, 34/255];     % Gray rgba(34, 34, 34, 1)
color_tolerance = [240/255, 228/255, 66/255]; % Yellow rgba(197, 166, 75, 1) - high visibility
color_timing = [0/255, 158/255, 115/255];     % Bluish green rgb(0, 158, 115) - colorblind safe

% Create data_save directory if it does not exist
if ~exist(data_save, 'dir')
    mkdir(data_save)
end
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir)
end
current_participant_name = 'test'; 
sessionPhaseNumber = 3; % default

% load data
subject_data = module_read_neurone(...
        fullfile(data_load, current_participant_name), ...
        sessionPhaseNumber = sessionPhaseNumber);


channel_labels = fieldnames(subject_data.signal)';
channel_labels(end-1:end) = [];
rs_EEG = [];
for i = 1:length(channel_labels)
    rs_EEG(:,i) = subject_data.signal.(channel_labels{i}).data;
end

%% Create EEGlab-compatible structure
EEG_data = eeg_emptyset();
EEG_data.data = rs_EEG';
EEG_data.srate = subject_data.properties.samplingRate;
EEG_data.nbchan = size(rs_EEG', 1);
EEG_data.pnts = size(rs_EEG', 2);
EEG_data.trials = 1;
EEG_data.xmin = 0;
EEG_data.xmax = (EEG_data.pnts - 1) / EEG_data.srate;
EEG_data.times = linspace(EEG_data.xmin, EEG_data.xmax, EEG_data.pnts) * 1000; % in ms

% Create chanlocs structure compatible with EEGlab
for i = 1:length(channel_labels)
    EEG_data.chanlocs(i).labels = channel_labels{i};
end

% Save data
save([data_save 'LAVA_EEG_rs_' current_participant_name '.mat'], "EEG_data", '-v7.3')

%% Collect all markers
% Extract all marker information from subject_data
triggers = subject_data.markers;

random_phases = readtable(data_csv);

BOSS_types = {'10', '20', '30'}; % BOSS trigger types 
NS_type = 'Stimulation';
NS_port = 'B';

BOSS_Triggers.types = triggers.type(ismember(triggers.type, BOSS_types));
BOSS_Triggers.times = triggers.time(ismember(triggers.type, BOSS_types));
BOSS_Triggers.times = round(BOSS_Triggers.times*1000);
BOSS_Triggers.phases = zeros(size(BOSS_Triggers.types));
NS_Triggers.times = triggers.time(ismember(triggers.port, NS_port) & ismember(triggers.type, NS_type));
NS_Triggers.times = round(NS_Triggers.times*1000);
NS_Triggers.phases = zeros(size(NS_Triggers.times));

random_phase_idx = 1; % Counter for type '30' triggers
for i = 1:length(BOSS_Triggers.types)
    switch BOSS_Triggers.types{i}
        case '10'
            BOSS_Triggers.phases(i) = pi;
        case '20'
            BOSS_Triggers.phases(i) = 0;
        case '30'
            BOSS_Triggers.phases(i) = random_phases.Var1(random_phase_idx);
            random_phase_idx = random_phase_idx + 1;
    end
end


% use the same trigger phases
if any(BOSS_Triggers.phases) && length(BOSS_Triggers.times) == length(NS_Triggers.times)
    NS_Triggers.phases = BOSS_Triggers.phases;
    fprintf("BOSS phases found or matching array sizes, copying over.\n")
else
    fprintf("No BOSS phases found, indexing over times.\n")

    random_phase_idx = 1;
    for i = 1:length(NS_Triggers.times)
        if 0 < i && i < 11
            NS_Triggers.phases(i) = pi;
    
        elseif 10 < i && i < 21
    
        else
            NS_Triggers.phases(i) = random_phases.phase(random_phase_idx);
            random_phase_idx = random_phase_idx + 1;
        end
    end
end

%% Preprocessing
hjorth.channel = 'C3';

switch hjorth.channel  % get surrounding variables for hjorth filter
    case 'C3'
        hjorth.electrodes = {'FC1', 'CP1', 'FC5', 'CP5'};
    case 'C4'
        hjorth.electrodes = {'FC2', 'CP2', 'FC6', 'CP6'};
end

down_fs = 500; % downsample frequency in Hz

EEG = prep_ds_dt_lap(EEG_data, hjorth, down_fs);


%% Bandpass filter

% Design the bandpass FIR filter
settings.filter = load("../data/filter_coeffs.mat").coeffs;

% Filter
data = filtfilt(settings.filter, 1, EEG.data')';

%% Phase estimation

% hilbert transform
hilb = hilbert(data);  % Applied to full data
phase.ground_truth = angle(hilb);

% get the indexes of trigger times in downsampled data
% Find nearest time points since downsampling creates discrete time points
NS_times_ms = NS_Triggers.times; % already in ms from line 55
BOSS_times_ms = BOSS_Triggers.times;
NS_inds_ds = zeros(size(NS_times_ms));
BOSS_inds_ds = zeros(size(BOSS_times_ms));

for i = 1:length(NS_times_ms)
    [~, NS_inds_ds(i)] = min(abs(EEG.times - NS_times_ms(i)));
end
fprintf('Found %d NS trigger indices in downsampled EEG data.\n', length(NS_inds_ds));
fprintf('Trigger time differences (ms):\n');
for i = 1:length(NS_inds_ds)  % Show first 5 as example
    fprintf(' BOSS Trigger %d: requested %.2f ms, found %.2f ms (diff: %.2f ms)\n', ...
        i, NS_times_ms(i), EEG.times(NS_inds_ds(i)), ...
        NS_times_ms(i) - EEG.times(NS_inds_ds(i)));
end


for i = 1:length(BOSS_times_ms)
    [~, BOSS_inds_ds(i)] = min(abs(EEG.times - BOSS_times_ms(i)));
end
fprintf('Found %d BOSS trigger indices in downsampled EEG data.\n', length(BOSS_inds_ds));
fprintf('Trigger time differences (ms):\n');
for i = 1:length(BOSS_inds_ds)  % Show first 5 as example
    fprintf('  BOSS Trigger %d: requested %.2f ms, found %.2f ms (diff: %.2f ms)\n', ...
        i, BOSS_times_ms(i), EEG.times(BOSS_inds_ds(i)), ...
        BOSS_times_ms(i) - EEG.times(BOSS_inds_ds(i)));
end


% check angle at trigger times 
phase.NS_at_triggers = phase.ground_truth(NS_inds_ds);
%phase.NS_at_triggers = mod(phase.NS_at_triggers + 2*pi, 2*pi); % convert to [0, 2pi)
phase.target = NS_Triggers.phases;
NS_phase_error = phase.NS_at_triggers - phase.target';
NS_mod_error = mod(NS_phase_error + pi, 2*pi) - pi; % wrap to [-pi, pi]
phase.BOSS_at_triggers = phase.ground_truth(BOSS_inds_ds);
%phase.BOSS_at_triggers = mod(phase.BOSS_at_triggers + 2*pi, 2*pi); % convert to [0, 2pi)
BOSS_phase_error = phase.BOSS_at_triggers - phase.target';
BOSS_mod_error = mod(BOSS_phase_error + pi, 2*pi) - pi; % wrap to [-pi, pi]

corr_BOSS_NS = corrcoef(BOSS_mod_error, NS_mod_error);
disp(["Correlation between NS and BOSS", corr_BOSS_NS(1,2)]);

NS_times_sec = NS_Triggers.times / 1000; % convert ms to seconds
BOSS_times_sec = BOSS_Triggers.times / 1000;

% get difference in timing per trial
diff_trigger = BOSS_times_sec - NS_times_sec;


% Group indices by condition
idx = 1:numel(NS_mod_error);

[G, group_labels] = findgroups(BOSS_Triggers.types);

grouped_NS = splitapply(@(i) {NS_mod_error(i)}, idx', G);
grouped_BOSS = splitapply(@(i) {BOSS_mod_error(i)}, idx', G);

% Compute means directly from numeric vectors
mean_NS_mod = cellfun(@(x) circ_mean(abs(x), [], 2), grouped_NS);
mean_BOSS_mod = cellfun(@(x) circ_mean(abs(x), [], 2), grouped_BOSS);

%% Statistical analysis

% Set RNG seed for reproducibility
rng(42);

% create x sampling grid
theta_grid = linspace(0, 2*pi, 720);

perms = 3000;

angles_NS = phase.NS_at_triggers;  % radians (0 tp 2pi)
angles_BOSS = phase.BOSS_at_triggers;

NS_pi_angles = angles_NS(G==1)';
NS_zero_angles = angles_NS(G==2)';
NS_random_angles = angles_NS(G==3)';
BOSS_pi_angles = angles_BOSS(G==1)';
BOSS_zero_angles = angles_BOSS(G==2)';
BOSS_random_angles = angles_BOSS(G==3)';


[test_NS.pi.params(1), test_NS.pi.params(2)] = circ_vmpar(NS_pi_angles);
[test_NS.zero.params(1), test_NS.zero.params(2)] = circ_vmpar(NS_zero_angles);
[test_NS.random.params(1), test_NS.random.params(2)] = circ_vmpar(NS_random_angles);
[test_BOSS.pi.params(1), test_BOSS.pi.params(2)] = circ_vmpar(BOSS_pi_angles);
[test_BOSS.zero.params(1), test_BOSS.zero.params(2)] = circ_vmpar(BOSS_zero_angles);
[test_BOSS.random.params(1), test_BOSS.random.params(2)] = circ_vmpar(BOSS_random_angles);

fprintf('\nFitted von Mises parameters:\n');
fprintf('   NS Target pi: mu = %.4f, kappa = %.4f\n', test_NS.pi.params(1), test_NS.pi.params(2));
fprintf('   NS Target 0: mu = %.4f, kappa = %.4f\n', test_NS.zero.params(1), test_NS.zero.params(2));
fprintf('   NS Target Random: mu = %.4f, kappa = %.4f\n', test_NS.random.params(1), test_NS.random.params(2));
fprintf('   BOSS Target pi: mu = %.4f, kappa = %.4f\n', test_BOSS.pi.params(1), test_BOSS.pi.params(2));
fprintf('   BOSS Target 0: mu = %.4f, kappa = %.4f\n', test_BOSS.zero.params(1), test_BOSS.zero.params(2));
fprintf('   BOSS Target Random: mu = %.4f, kappa = %.4f\n', test_BOSS.random.params(1), test_BOSS.random.params(2));

% 1. test whehter the agnular data follows von Mises distribution

test_NS.pi.theor_dist = circ_vmrnd(test_NS.pi.params(1), test_NS.pi.params(2), length(NS_pi_angles));
test_NS.zero.theor_dist = circ_vmrnd(test_NS.zero.params(1), test_NS.zero.params(2), length(NS_zero_angles));
test_NS.random.theor_dist = circ_vmrnd(test_NS.random.params(1), test_NS.random.params(2), length(NS_random_angles));
test_BOSS.pi.theor_dist = circ_vmrnd(test_BOSS.pi.params(1), test_BOSS.pi.params(2), length(BOSS_pi_angles));
test_BOSS.zero.theor_dist = circ_vmrnd(test_BOSS.zero.params(1), test_BOSS.zero.params(2), length(BOSS_zero_angles));
test_BOSS.random.theor_dist = circ_vmrnd(test_BOSS.random.params(1), test_BOSS.random.params(2), length(BOSS_random_angles));

% theoretical data for KDE comparison
kappa_pi = max(1,(test_NS.pi.params(2) + test_BOSS.pi.params(2)) / 2);      % Average kappa for pi target
kappa_zero = max(1,(test_NS.zero.params(2) + test_BOSS.zero.params(2)) / 2); % Average kappa for zero target
kappa_random = 0;                     % Keep as 0 for uniform (random) distribution

fprintf('\n(Averaged) kappa values for theoretical distributions:\n');
fprintf('   Target pi: kappa = %.4f\n', kappa_pi);
fprintf('   Target 0: kappa = %.4f\n', kappa_zero);
fprintf('   Target Random: kappa = %.4f\n', kappa_random);



test_pi.theor_data = circ_vmpdf(theta_grid, pi, kappa_pi)';       % Use averaged kappa for pi
test_zero.theor_data = circ_vmpdf(theta_grid, 0, kappa_zero)';    % Use averaged kappa for zero
test_random.theor_data = circ_vmpdf(theta_grid, 0, kappa_random)'; % Kappa=0 for uniform


% test whether angles_NS follows von Mises distribution
[test_NS.pi.p_val.wu2_perm, ~] = watsons_U2_perm_test(NS_pi_angles, test_NS.pi.theor_dist, perms);
[test_NS.zero.p_val.wu2_perm, ~] = watsons_U2_perm_test(NS_zero_angles, test_NS.zero.theor_dist, perms);
[test_NS.random.p_val.wu2_perm, ~] = watsons_U2_perm_test(NS_random_angles, test_NS.random.theor_dist, perms);
% test whether angles_BOSS follows von Mises distribution
[test_BOSS.pi.p_val.wu2_perm, ~] = watsons_U2_perm_test(BOSS_pi_angles, test_BOSS.pi.theor_dist, perms);
[test_BOSS.zero.p_val.wu2_perm, ~] = watsons_U2_perm_test(BOSS_zero_angles, test_BOSS.zero.theor_dist, perms);
[test_BOSS.random.p_val.wu2_perm, ~] = watsons_U2_perm_test(BOSS_random_angles, test_BOSS.random.theor_dist, perms);

fprintf('\nCircular statistics results:\n');
fprintf('   NS Phase Fit to von Mises: \n');
fprintf('       Target π: p = %.4f, mu = %.4f, kappa = %.4f\n', test_NS.pi.p_val.wu2_perm, test_NS.pi.params(1), test_NS.pi.params(2));
fprintf('       Target 0: p = %.4f, mu = %.4f, kappa = %.4f\n', test_NS.zero.p_val.wu2_perm, test_NS.zero.params(1), test_NS.zero.params(2));
fprintf('       Target Random: p = %.4f, mu = %.4f, kappa = %.4f\n', test_NS.random.p_val.wu2_perm, test_NS.random.params(1), test_NS.random.params(2));
fprintf('   BOSS Phase Fit to von Mises: \n');
fprintf('       Target π: p = %.4f, mu = %.4f, kappa = %.4f\n', test_BOSS.pi.p_val.wu2_perm, test_BOSS.pi.params(1), test_BOSS.pi.params(2));
fprintf('       Target 0: p = %.4f, mu = %.4f, kappa = %.4f\n', test_BOSS.zero.p_val.wu2_perm, test_BOSS.zero.params(1), test_BOSS.zero.params(2));
fprintf('       Target Random: p = %.4f, mu = %.4f, kappa = %.4f\n', test_BOSS.random.p_val.wu2_perm, test_BOSS.random.params(1), test_BOSS.random.params(2));

% Test 3: Check whether phase distributions have a preferred direction
[test_NS.pi.p_val.rayleigh, ~] = circ_rtest(NS_pi_angles);
[test_NS.zero.p_val.rayleigh, ~] = circ_rtest(NS_zero_angles);
[test_BOSS.pi.p_val.rayleigh, ~] = circ_rtest(BOSS_pi_angles);
[test_BOSS.zero.p_val.rayleigh, ~] = circ_rtest(BOSS_zero_angles);
[test_NS.random.p_val.rayleigh, ~] = circ_rtest(NS_random_angles);
[test_BOSS.random.p_val.rayleigh, ~] = circ_rtest(BOSS_random_angles);
fprintf('\nRayleigh test for non-uniformity around target phase:\n');
fprintf('   NS Target pi: p = %.4f\n', test_NS.pi.p_val.rayleigh);
fprintf('   NS Target 0: p = %.4f\n', test_NS.zero.p_val.rayleigh);
fprintf('   BOSS Target pi: p = %.4f\n', test_BOSS.pi.p_val.rayleigh);
fprintf('   BOSS Target 0: p = %.4f\n', test_BOSS.zero.p_val.rayleigh);
fprintf('   NS Target Random: p = %.4f\n', test_NS.random.p_val.rayleigh);
fprintf('   BOSS Target Random: p = %.4f\n', test_BOSS.random.p_val.rayleigh);

% Test 4: Test for a target mean of our empirical data
[test_NS.pi.p_val.v_test, ~] = circ_vtest(NS_pi_angles, pi);
[test_NS.zero.p_val.v_test, ~] = circ_vtest(NS_zero_angles, 0);
[test_BOSS.pi.p_val.v_test, ~] = circ_vtest(BOSS_pi_angles, pi);
[test_BOSS.zero.p_val.v_test, ~] = circ_vtest(BOSS_zero_angles, 0);
fprintf('\nV test for mean phase equal to target phase:\n');
fprintf('   NS Target pi: p = %.4f\n', test_NS.pi.p_val.v_test);
fprintf('   NS Target 0: p = %.4f\n', test_NS.zero.p_val.v_test);
fprintf('   BOSS Target pi: p = %.4f\n', test_BOSS.pi.p_val.v_test);
fprintf('   BOSS Target 0: p = %.4f\n', test_BOSS.zero.p_val.v_test);

% FOR CONVENIENCE!!! TRANSFER THE P VALUE OF THE RAYLEIGH TEST FOR THE RANDOM PHASE TARGET TO V TEST
test_NS.random.p_val.v_test = test_NS.random.p_val.rayleigh;
test_BOSS.random.p_val.v_test = test_BOSS.random.p_val.rayleigh;


% Test 5: Check whether each distributed phase group differs between NS and BOSS
% Using Watson's U2 test with permutation
[test_pi.BOSS_NS.p_val, ~] = watsons_U2_perm_test(NS_pi_angles, BOSS_pi_angles, perms);
[test_zero.BOSS_NS.p_val, ~] = watsons_U2_perm_test(NS_zero_angles, BOSS_zero_angles, perms);
[test_random.BOSS_NS.p_val, ~] = watsons_U2_perm_test(NS_random_angles, BOSS_random_angles, perms);

fprintf('\nWatson U2 test between NS and BOSS:\n');
fprintf('   Target pi: p = %.4f\n', test_pi.BOSS_NS.p_val);
fprintf('   Target 0: p = %.4f\n', test_zero.BOSS_NS.p_val);
fprintf('   Target Random: p = %.4f\n', test_random.BOSS_NS.p_val);

% get correlation between NS and BOSS Data: is this a good quanitiy measure of similarity?
[test_pi.ksd_corr_NS_BOSS(1), test_pi.ksd_corr_NS_BOSS(2)] = circ_corrcc(NS_pi_angles, BOSS_pi_angles);
[test_zero.ksd_corr_NS_BOSS(1), test_zero.ksd_corr_NS_BOSS(2)] = circ_corrcc(NS_zero_angles, BOSS_zero_angles);
[test_random.ksd_corr_NS_BOSS(1), test_random.ksd_corr_NS_BOSS(2)] = circ_corrcc(NS_random_angles, BOSS_random_angles);
fprintf('\nKDE correlation between NS and BOSS:\n');
fprintf('   Target pi: r = %.4f, with p = %.6f\n', test_pi.ksd_corr_NS_BOSS(1), test_pi.ksd_corr_NS_BOSS(2));
fprintf('   Target 0: r = %.4f, with p = %.6f\n', test_zero.ksd_corr_NS_BOSS(1), test_zero.ksd_corr_NS_BOSS(2));
fprintf('   Target Random: r = %.4f, with p = %.6f\n', test_random.ksd_corr_NS_BOSS(1), test_random.ksd_corr_NS_BOSS(2));

% get conf_interval for mean phase
% test_NS.pi.conf_mean = circ_confmean(NS_pi_angles');
% test_NS.zero.conf_mean = circ_confmean(NS_zero_angles');
% test_BOSS.pi.conf_mean = circ_confmean(BOSS_pi_angles');
% test_BOSS.zero.conf_mean = circ_confmean(BOSS_zero_angles');
test_NS.pi.conf_mean = [pi/8, 3*pi/8]; % placeholder
test_NS.zero.conf_mean = [-pi/8, pi/8]; % placeholder
test_BOSS.pi.conf_mean = [pi/6, 2*pi/5]; % placeholder
test_BOSS.zero.conf_mean = [-pi/6, pi/6]; % placeholder


fprintf('\nConfidence intervals for mean phase:\n');
fprintf('   NS Target pi: [%.4f, %.4f]\n', test_NS.pi.conf_mean(1), test_NS.pi.conf_mean(2));
fprintf('   NS Target 0: [%.4f, %.4f]\n', test_NS.zero.conf_mean(1), test_NS.zero.conf_mean(2));
fprintf('   BOSS Target pi: [%.4f, %.4f]\n', test_BOSS.pi.conf_mean(1), test_BOSS.pi.conf_mean(2));
fprintf('   BOSS Target 0: [%.4f, %.4f]\n', test_BOSS.zero.conf_mean(1), test_BOSS.zero.conf_mean(2));

% make a KDE per condition and algorithm
% get the theoretical
test_NS.pi.ksd_estimate = circ_ksdensity(NS_pi_angles, theta_grid);
test_NS.zero.ksd_estimate = circ_ksdensity(NS_zero_angles, theta_grid);
test_NS.random.ksd_estimate = circ_ksdensity(NS_random_angles, theta_grid);
test_BOSS.pi.ksd_estimate = circ_ksdensity(BOSS_pi_angles, theta_grid);
test_BOSS.zero.ksd_estimate = circ_ksdensity(BOSS_zero_angles, theta_grid);
test_BOSS.random.ksd_estimate = circ_ksdensity(BOSS_random_angles, theta_grid);

% get the similarity between empirical and theoretical distributions
[test_NS.pi.ksd_corr(1), test_NS.pi.ksd_corr(2)] = circ_corrcc(test_NS.pi.ksd_estimate, test_pi.theor_data);
[test_NS.zero.ksd_corr(1), test_NS.zero.ksd_corr(2)] = circ_corrcc(test_NS.zero.ksd_estimate, test_zero.theor_data);
[test_NS.random.ksd_corr(1), test_NS.random.ksd_corr(2)] = circ_corrcc(test_NS.random.ksd_estimate, test_random.theor_data);
[test_BOSS.pi.ksd_corr(1), test_BOSS.pi.ksd_corr(2)] = circ_corrcc(test_BOSS.pi.ksd_estimate, test_pi.theor_data);
[test_BOSS.zero.ksd_corr(1), test_BOSS.zero.ksd_corr(2)] = circ_corrcc(test_BOSS.zero.ksd_estimate, test_zero.theor_data);
[test_BOSS.random.ksd_corr(1), test_BOSS.random.ksd_corr(2)] = circ_corrcc(test_BOSS.random.ksd_estimate, test_random.theor_data);

fprintf('\nKDE correlation with theoretical distribution:\n');
fprintf('   NS Target pi: r = %.4f, with p = %.4f\n', test_NS.pi.ksd_corr(1), test_NS.pi.ksd_corr(2));
fprintf('   NS Target 0: r = %.4f, with p = %.4f\n', test_NS.zero.ksd_corr(1), test_NS.zero.ksd_corr(2));
fprintf('   NS Target Random: r = %.4f, with p = %.4f\n', test_NS.random.ksd_corr(1), test_NS.random.ksd_corr(2));
fprintf('   BOSS Target pi: r = %.4f, with p = %.4f\n', test_BOSS.pi.ksd_corr(1), test_BOSS.pi.ksd_corr(2));
fprintf('   BOSS Target 0: r = %.4f, with p = %.4f\n', test_BOSS.zero.ksd_corr(1), test_BOSS.zero.ksd_corr(2));
fprintf('   BOSS Target Random: r = %.4f, with p = %.4f\n', test_BOSS.random.ksd_corr(1), test_BOSS.random.ksd_corr(2));  

% compute the overlapping coefficient between empircial and theoretical distributions
test_NS.pi.ksd_ovl = trapz(theta_grid, min([test_NS.pi.ksd_estimate; test_pi.theor_data]));
test_NS.zero.ksd_ovl = trapz(theta_grid, min([test_NS.zero.ksd_estimate; test_zero.theor_data]));
test_NS.random.ksd_ovl = trapz(theta_grid, min([test_NS.random.ksd_estimate; test_random.theor_data]));
test_BOSS.pi.ksd_ovl = trapz(theta_grid, min([test_BOSS.pi.ksd_estimate; test_pi.theor_data]));
test_BOSS.zero.ksd_ovl = trapz(theta_grid, min([test_BOSS.zero.ksd_estimate; test_zero.theor_data]));
test_BOSS.random.ksd_ovl = trapz(theta_grid, min([test_BOSS.random.ksd_estimate; test_random.theor_data])); 
fprintf('\nKDE overlapping coefficient with theoretical distribution:\n');
fprintf('   NS Target pi: OVL = %.4f\n', test_NS.pi.ksd_ovl);
fprintf('   NS Target 0: OVL = %.4f\n', test_NS.zero.ksd_ovl);
fprintf('   NS Target Random: OVL = %.4f\n', test_NS.random.ksd_ovl);
fprintf('   BOSS Target pi: OVL = %.4f\n', test_BOSS.pi.ksd_ovl);
fprintf('   BOSS Target 0: OVL = %.4f\n', test_BOSS.zero.ksd_ovl);
fprintf('   BOSS Target Random: OVL = %.4f\n', test_BOSS.random.ksd_ovl);


%% Compute the SNR PSD of the EEG data
fs = down_fs; % Sampling frequency after downsampling
eeg_signal = EEG.data; % Use channel 1, or change as needed


% Compute Power Spectral Density (PSD) using Welch's method
window = hann(4*fs); % 4-second window
noverlap = length(window) - 10;
nfft = max(2^nextpow2(length(eeg_signal)), 1024);

[pxx, f] = pwelch(eeg_signal, window, noverlap, nfft, fs);

% SNR calculation: subtract 1/f fit from spectrum
mask_1f = (f >= 0.5 & f <= 8) | (f >= 30 & f <= 45);
mask_mu = f >= 8 & f <= 15;
mask_full = f >= 0.5 & f < 45;

% Fit 1/f noise (log-log fit)
fit_freqs = f(mask_1f);
fit_power = 10*log10(pxx(mask_1f));
p = polyfit(log10(fit_freqs), fit_power, 1);
yfit = polyval(p, log10(f(mask_full)));

% Subtract 1/f fit from spectrum
new_spec = 10*log10(pxx(mask_full)) - yfit;
freqs = f(mask_full);

[H, f_filter] = freqz(settings.filter, 1, length(freqs), fs); % get frequency response of filter

y_filter =(H.*conj(H))'; % get power response of filter

%% Plot 0: SNR over Frequency
snr_threshold = 5; % SNR threshold in dB


figure("Name", "SNR over Frequency", 'Units','pixels','Position',[400 400 3400 1800]);
hold on;
yline(snr_threshold, '-', 'DisplayName', 'SNR Threshold', 'LineWidth', 1, 'Color', 'r');
ylim([min(new_spec) - 1, max(max(new_spec), snr_threshold) + 5]);
ylabel('SNR [dB]');

plot(f_filter, y_filter*2, '--', 'LineWidth', 0.75, 'Color', '#4e4e4e', 'DisplayName', sprintf('Power Response of\nBandpass Filter (scaled)'));
plot(freqs, new_spec, 'LineWidth', 1.5, 'DisplayName', 'SNR (Spectrum minus 1/f)', 'Color', '#0072B2');

xlim([0 max(freqs)]);
xlabel('Frequency (Hz)');
title('SNR over Frequency', 'FontSize', 12);
legend;
grid on;
hold off;

%% Plot 1: Phase error and trigger times
figure('Units','pixels','Position',[200 200 2600 2000]);

% Define relative positions [left bottom width height]
top_pos = [0.05 0.55 0.9 0.4];    % top plot (full width, taller)
bottom_left_pos = [0.05 0.06 0.3 0.4];  % bottom-left plot
bottom_right_pos = [0.4 0.06 0.55 0.4]; % bottom-right plot

% top
ax1 = axes('Position', top_pos);
hold on;

% Plot scatter first to establish data range
scatter(idx, NS_mod_error, 15, color_NS, 'filled');
scatter(idx, BOSS_mod_error, 15, color_BOSS, 'filled');
ylim([-pi-0.1 pi+0.1]);


max_x = length(NS_mod_error) + 1;
% Add tolerance band (in background) spanning only data range
tolerance = pi/40;
patch([-1 max_x max_x -1], ...
     [-tolerance -tolerance tolerance tolerance], ...
     color_tolerance, 'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
     'DisplayName', 'Tolerance (±π/40)');
plot([-1 max_x], [tolerance tolerance], '-', 'Color', color_tolerance, 'HandleVisibility', 'off');
plot([-1 max_x], [-tolerance -tolerance], '-', 'Color', color_tolerance, 'HandleVisibility', 'off');

% Bring scatter to front
uistack(findobj(gca, 'Type', 'scatter'), 'top');

x_offset = -0.8; % offset to the left in seconds

% detect conditiion changes in BOSS_Trrigers.types
cond_changes = [1; find(diff(cellfun(@str2double, BOSS_Triggers.types)) ~= 0) + 1; length(BOSS_Triggers.types) + 1];


    % Determine condition name and color based on target value
condition_name = 'π (trough)';
xline(cond_changes(1) + x_offset, '-', condition_name, ...
    'LineWidth', 1, 'LabelVerticalAlignment', 'bottom', 'Color', color_line);

condition_name = '0 (peak)';
xline(cond_changes(2) + x_offset, '-', condition_name, ...
    'LineWidth', 1, 'LabelVerticalAlignment', 'bottom', 'Color', color_line);



condition_name = 'Random';
xline(cond_changes(3) + x_offset, '-', condition_name, ...
    'LineWidth', 1, 'LabelVerticalAlignment', 'bottom', 'Color', color_line);

xlim([-0.5, max_x]);
xlabel('Trial Number');
ylabel('Phase Error (radians)');
legend('Phase Tolerance', 'Phase Error NS', 'Phase Error NS', 'Location', 'southeast');
title('Phase Error at Trigger Times');      
grid on;
hold off;

% Bottom-left: Mean Phase Error
ax2 = axes('Position', bottom_left_pos);
b = bar(ax2, [mean_NS_mod, mean_BOSS_mod]);
b(1).FaceColor = color_NS;
b(2).FaceColor = color_BOSS;
x = {"π", "0", "Random"};
set(ax2, 'XTick', 1:numel(x), 'XTickLabel', x);
xlabel(ax2,"Phase Target");
ylabel(ax2,"Mean Error (rad)");
title(ax2,"Mean Phase Error (absolute)");
legend(ax2,'NS','BOSS','Location','southeast');

% Bottom-right: Trigger differences
ax3 = axes('Position', bottom_right_pos);
plot(ax3, 1:numel(diff_trigger), diff_trigger, '-o', 'Color', color_timing, 'MarkerFaceColor', color_timing);
xlabel(ax3,"Trial Number");
xlim([0, length(diff_trigger) + 1])
ylabel(ax3,"Trigger difference (s)");
title(ax3,"Timing differences (BOSS - NS)");
grid(ax3,'on');

saveas(gcf, fullfile(fig_dir, ['phaseerror_and_triggertimes_' current_participant_name '.png']));

%% Plot 2: 2x3 polar-grid (rows: NS, BOSS) × (cols: π, 0, Random)
% Prepare target mapping and angle arrays
targets = {'π','0','Random'};
target_fields = {'pi', 'zero', 'random'};
target_angles = [pi, 0, NaN]; % NaN for Random (no target line)
alg_labels = {'NS','BOSS'};
% Determine global max for radial limit
all_angles = [NS_pi_angles; NS_zero_angles; NS_random_angles; ...
              BOSS_pi_angles; BOSS_zero_angles; BOSS_random_angles];
[counts, edges] = histcounts(all_angles, 30, 'Normalization', 'probability');
global_max = max(counts);

% Create figure and layout
figure('Units','pixels','Position',[200 200 3600 2000]);
left_margin = 0.06; right_margin = 0.01; top_margin = 0.01; bottom_margin = 0.05;
cols = 3; rows = 2;
width = (1.05 - left_margin - right_margin) / cols;
height = (1.02 - top_margin - bottom_margin) / rows;
fontsize = 10;
for col = 1:cols
    for row = 1:rows
        % compute normalized position: origin bottom-left
        left = left_margin + (col-1) * width;
        % rows: row=1 should be top (NS), row=2 bottom (BOSS)
        vertical_gap = -0.035; % adjust gap between rows
        height_adj = height * 0.8; 
        bottom = 1 - top_margin - row * height - (row-1)*vertical_gap;
        ax = polaraxes('Position', [left, bottom, width*0.8, height_adj]);


        if row == 1
            data = angles_NS(G==col);
            mu_hat = test_NS.(target_fields{col}).params(1);
            col_color = color_NS;
            p_val = test_NS.(target_fields{col}).p_val.v_test;
        else
            data = angles_BOSS(G==col);
            mu_hat = test_BOSS.(target_fields{col}).params(1);
            col_color = color_BOSS;
            p_val = test_BOSS.(target_fields{col}).p_val.v_test;
        end

        
        h = polarhistogram(ax, data, 30, 'Normalization', 'probability', ...
            'FaceColor', col_color, 'FaceAlpha', 0.75, 'EdgeColor', 'none', Facealpha = 0.7);

        % set radial limit consistent across panels
        rlim(ax, [0 global_max*1.1]);
        thetaticks(ax, 0:45:315);
        thetaticklabels(ax, {'0°','45°','90°','135°','180°','225°','270°','315°'});
        ax.ThetaDir = 'clockwise';
        ax.ThetaZeroLocation = 'right';

        % draw a radial line for target if not Random
        if ~isnan(target_angles(col))
            t = target_angles(col);
            t2 = mod(t, 2*pi);
            hold(ax, 'on');
            polarplot(ax, [t2 t2], [0.003 global_max*1.1], '--k', 'LineWidth', 1.0);
            
            hold(ax, 'off');
            % confidence interval shading
            if row == 1
                conf_int = test_NS.(target_fields{col}).conf_mean;
                color_conf = color_NS*0.7;
            else
                conf_int = test_BOSS.(target_fields{col}).conf_mean;
                color_conf = color_BOSS*0.6;
            end
            % Define the angular range between the two confidence bounds
            theta_arc = linspace(conf_int(1), conf_int(2), 100); % angles for the arc
            r_arc = global_max*1.1;     % fixed radius at outer edge

            hold(ax, 'on');
            polarplot(ax, theta_arc, r_arc * ones(size(theta_arc)), ...
                'Color', color_conf, 'LineWidth', 4);
            polarplot(ax, [mu_hat mu_hat], [0.003 global_max*1.1], '-', 'Color', color_conf, 'LineWidth', 1.25);
            
            
            hold(ax, 'off');  
        end
        
            text(ax, -pi, global_max*1.5, sprintf('%s — %s\n^{(p = %.3f)}', alg_labels{row}, targets{col}, p_val), ...
            'Rotation', 90, ...         % Rotate text vertically
            'HorizontalAlignment', 'center', ...
            'FontSize', fontsize, 'FontWeight', 'bold');     
    end

    % add the circular corr coeff between NS and BOSS for each target
    circ_correff = eval(['test_' target_fields{col} '.ksd_corr_NS_BOSS(1)']);
    % just below the 90° — render 'circ' as a lowercase subscript
    text(ax, pi/2, global_max*1.35, sprintf('r_{%s} = %.2f', 'circ', circ_correff), ...
        'HorizontalAlignment', 'center', ...
        'FontSize', fontsize, 'FontWeight', 'bold', ...
        'Interpreter', 'tex');


end


sgtitle('Phase distributions by target and algorithm', ...
    'FontSize', fontsize + 4, 'FontWeight', 'bold');

% Save figure
saveas(gcf, fullfile(fig_dir, ['phase_distributions_polar_' current_participant_name '.png']));

%% Plot 3: KDE estimates of phase distributions
figure('Units','pixels','Position',[400 400 4500 1700]);

line_width = 1.1;
alpha_hist = 0.2;

% polar plot 1 for pi target
ax1 = polaraxes('Position', [0.05, 0.05, 0.25, 0.85]);
hold(ax1, 'on');
h_pi_NS = polarhistogram(ax1, NS_pi_angles, 30, 'Normalization', 'pdf', ...
    'FaceColor', color_NS, 'FaceAlpha', alpha_hist, 'EdgeColor', 'none', 'HandleVisibility', 'off');
h_pi_BOSS = polarhistogram(ax1, BOSS_pi_angles, 30, 'Normalization', 'pdf', ...
    'FaceColor', color_BOSS, 'FaceAlpha', alpha_hist, 'EdgeColor', 'none','HandleVisibility', 'off');
kde_scale_1 = max([h_pi_NS.Values, h_pi_BOSS.Values]) / max([test_NS.pi.ksd_estimate, test_BOSS.pi.ksd_estimate, test_pi.theor_data]);  
polarplot(ax1, theta_grid, test_NS.pi.ksd_estimate*kde_scale_1, 'LineWidth', line_width, 'Color', color_NS);
polarplot(ax1, theta_grid, test_BOSS.pi.ksd_estimate*kde_scale_1, 'LineWidth', line_width, 'Color', color_BOSS);
polarplot(ax1, theta_grid, test_pi.theor_data*kde_scale_1, 'LineWidth', line_width, 'Color', color_line);
hold(ax1, 'off');
pi_max = max([h_pi_NS.Values, h_pi_BOSS.Values, max(test_pi.theor_data*kde_scale_1)]);
text(ax1, -pi/2, pi_max*1.4+0.09, sprintf('Target π\nr_{circ, NS} = %.3f / r_{circ, BOSS} = %.3f', test_NS.pi.ksd_corr(1), test_BOSS.pi.ksd_corr(1)), ...
            'HorizontalAlignment', 'center', ...
            'FontSize', fontsize, 'FontWeight', 'bold');   

% polar plot 2 for zero target
ax2 = polaraxes('Position', [0.38, 0.05, 0.25, 0.85]);
hold(ax2, 'on');
h_zero_NS = polarhistogram(ax2, NS_zero_angles, 30, 'Normalization', 'pdf', ...
    'FaceColor', color_NS, 'FaceAlpha', alpha_hist, 'EdgeColor', 'none','HandleVisibility', 'off');
h_zero_BOSS = polarhistogram(ax2, BOSS_zero_angles, 30, 'Normalization', 'pdf', ...
    'FaceColor', color_BOSS, 'FaceAlpha', alpha_hist, 'EdgeColor', 'none','HandleVisibility', 'off');
kde_scale_2 = max([h_zero_NS.Values, h_zero_BOSS.Values]) / max([test_NS.zero.ksd_estimate, test_BOSS.zero.ksd_estimate, test_zero.theor_data]);  
polarplot(ax2, theta_grid, test_NS.zero.ksd_estimate*kde_scale_2, 'LineWidth', line_width, 'Color', color_NS);
polarplot(ax2, theta_grid, test_BOSS.zero.ksd_estimate*kde_scale_2, 'LineWidth', line_width, 'Color', color_BOSS);
polarplot(ax2, theta_grid, test_zero.theor_data*kde_scale_2, 'LineWidth', line_width, 'Color', color_line); 
hold(ax2, 'off');
zero_max = max([h_zero_NS.Values, h_zero_BOSS.Values, max(test_zero.theor_data*kde_scale_2)]);
text(ax2, -pi/2, zero_max*1.4, sprintf('Target 0\nr_{circ, NS} = %.3f / r_{circ, BOSS} = %.3f', test_NS.zero.ksd_corr(1), test_BOSS.zero.ksd_corr(1)), ...
            'HorizontalAlignment', 'center', ...
            'FontSize', fontsize, 'FontWeight', 'bold');

% polar plot 3 for random target
ax3 = polaraxes('Position', [0.71, 0.05, 0.25, 0.85]);
hold(ax3, 'on');
h_random_NS = polarhistogram(ax3, NS_random_angles, 30, 'Normalization', 'pdf', ...
    'FaceColor', color_NS, 'FaceAlpha', alpha_hist, 'EdgeColor', 'none','HandleVisibility', 'off');
h_random_BOSS = polarhistogram(ax3, BOSS_random_angles, 30, 'Normalization', 'pdf', ...
    'FaceColor', color_BOSS, 'FaceAlpha', alpha_hist, 'EdgeColor', 'none','HandleVisibility', 'off');
kde_scale_3 = max([h_random_NS.Values, h_random_BOSS.Values]) / max([test_NS.random.ksd_estimate, test_BOSS.random.ksd_estimate, test_random.theor_data]);  
polarplot(ax3, theta_grid, test_NS.random.ksd_estimate*kde_scale_3, 'LineWidth', line_width, 'Color', color_NS);
polarplot(ax3, theta_grid, test_BOSS.random.ksd_estimate*kde_scale_3, 'LineWidth', line_width, 'Color', color_BOSS);
polarplot(ax3, theta_grid, test_random.theor_data*kde_scale_3, 'LineWidth', line_width, 'Color', color_line); 
hold(ax3, 'off');
random_max = max([h_random_NS.Values, h_random_BOSS.Values, max(test_random.theor_data*kde_scale_3)]);
text(ax3, -pi/2, random_max*1.4, sprintf('Target Random\nr_{circ, NS} = %.3f / r_{circ, BOSS} = %.3f', test_NS.random.ksd_corr(1), test_BOSS.random.ksd_corr(1)), ...
            'HorizontalAlignment', 'center', ...
            'FontSize', fontsize, 'FontWeight', 'bold');

sgtitle('KDE of Phase Distributions at Trigger Times', 'FontSize', fontsize + 4, 'FontWeight', 'bold');
legend({'NS', 'BOSS', 'Theoretical'}, 'Location', [0.46, 0.85, 0.1, 0.03], 'Orientation', 'horizontal');
legend boxoff;

% Save figure
saveas(gcf, fullfile(fig_dir, ['phase_distributions_kde_' current_participant_name '.png']));

%% Export statistical results

% export mat
stat_results = struct();
stat_results.test_NS = test_NS;
stat_results.test_BOSS = test_BOSS;
stat_results.test_pi = test_pi;
stat_results.test_zero = test_zero;
stat_results.test_random = test_random;
save(fullfile(fig_dir, ['circular_statistics_results_' current_participant_name '.mat']), 'stat_results');
fprintf('Circular statistics results saved to %s\n', fullfile(fig_dir, ['circular_statistics_results_' current_participant_name '.mat']));

% export csv
results_table = table({'NS Target π'; 'NS Target 0'; 'NS Target Random'; ...
                       'BOSS Target π'; 'BOSS Target 0'; 'BOSS Target Random'}, ...
                      [test_NS.pi.params(1); test_NS.zero.params(1); test_NS.random.params(1); ...
                       test_BOSS.pi.params(1); test_BOSS.zero.params(1); test_BOSS.random.params(1)], ...
                      [test_NS.pi.params(2); test_NS.zero.params(2); test_NS.random.params(2); ...
                       test_BOSS.pi.params(2); test_BOSS.zero.params(2); test_BOSS.random.params(2)], ...
                      [test_NS.pi.p_val.wu2_perm; test_NS.zero.p_val.wu2_perm; test_NS.random.p_val.wu2_perm; ...
                       test_BOSS.pi.p_val.wu2_perm; test_BOSS.zero.p_val.wu2_perm; test_BOSS.random.p_val.wu2_perm], ...
                      [test_NS.pi.p_val.rayleigh; test_NS.zero.p_val.rayleigh; test_NS.random.p_val.rayleigh; ...
                       test_BOSS.pi.p_val.rayleigh; test_BOSS.zero.p_val.rayleigh; test_BOSS.random.p_val.rayleigh], ...
                      [test_NS.pi.p_val.v_test; test_NS.zero.p_val.v_test; NaN; ...
                       test_BOSS.pi.p_val.v_test; test_BOSS.zero.p_val.v_test; NaN], ...
                      [test_pi.BOSS_NS.p_val; test_zero.BOSS_NS.p_val; test_random.BOSS_NS.p_val; ...
                       NaN; NaN; NaN], ...
                      [test_NS.pi.ksd_corr(1); test_NS.zero.ksd_corr(1); test_NS.random.ksd_corr(1); ...
                       test_BOSS.pi.ksd_corr(1); test_BOSS.zero.ksd_corr(1); test_BOSS.random.ksd_corr(1)], ...
                      [test_NS.pi.ksd_corr(2); test_NS.zero.ksd_corr(2); test_NS.random.ksd_corr(2); ...
                       test_BOSS.pi.ksd_corr(2); test_BOSS.zero.ksd_corr(2); test_BOSS.random.ksd_corr(2)], ...
                       [test_NS.pi.ksd_ovl; test_NS.zero.ksd_ovl; test_NS.random.ksd_ovl; ...
                       test_BOSS.pi.ksd_ovl; test_BOSS.zero.ksd_ovl; test_BOSS.random.ksd_ovl], ...
                      [test_pi.ksd_corr_NS_BOSS(1); test_zero.ksd_corr_NS_BOSS(1); test_random.ksd_corr_NS_BOSS(1); ...
                       NaN; NaN; NaN], ...
                      [test_pi.ksd_corr_NS_BOSS(2); test_zero.ksd_corr_NS_BOSS(2); test_random.ksd_corr_NS_BOSS(2); ...
                       NaN; NaN; NaN], ...
                      'VariableNames', {'Condition', 'Mu', 'Kappa', 'p_value_Watson_U2', 'p_value_Rayleigh', 'p_value_V_test', ...
                                        'p_value_Watson_U2_NS_vs_BOSS','OVC_with_theor' 'KDE_corr_theor_r', 'KDE_corr_theor_p', ...
                                        'KDE_corr_NS_BOSS_r', 'KDE_corr_NS_BOSS_p'});

writetable(results_table, fullfile(fig_dir, ['circular_statistics_results_' current_participant_name '.csv']));
fprintf('Circular statistics results exported to %s\n', fullfile(fig_dir, ['circular_statistics_results_' current_participant_name '.csv']));

%% Generate LaTeX table from the results
latex_str = '\\begin{tabular}{|l|r|r|r|r|r|r|r|r|r|r|}\n\\hline\n';
latex_str = [latex_str 'Condition & Mu & Kappa & p Watson U2 & p Rayleigh & p V test & p Watson U2 NS vs BOSS & KDE corr theor r & KDE corr theor p & KDE corr NS BOSS r & KDE corr NS BOSS p \\\\\n\\hline\n'];

for i = 1:size(results_table, 1)
    row = results_table(i,:);
    latex_str = [latex_str row.Condition{1} ' & ' sprintf('%.4f', row.Mu) ' & ' sprintf('%.4f', row.Kappa) ' & ' sprintf('%.4f', row.p_value_Watson_U2) ' & '];
    if isnan(row.p_value_Rayleigh)
        latex_str = [latex_str 'N/A & '];
    else
        latex_str = [latex_str sprintf('%.4f', row.p_value_Rayleigh) ' & '];
    end
    if isnan(row.p_value_V_test)
        latex_str = [latex_str 'N/A & '];
    else
        latex_str = [latex_str sprintf('%.4f', row.p_value_V_test) ' & '];
    end
    latex_str = [latex_str sprintf('%.4f', row.p_value_Watson_U2_NS_vs_BOSS) ' & ' sprintf('%.4f', row.KDE_corr_theor_r) ' & ' sprintf('%.6f', row.KDE_corr_theor_p) ' & '];
    if isnan(row.KDE_corr_NS_BOSS_r)
        latex_str = [latex_str 'N/A & N/A \\\\\n'];
    else
        latex_str = [latex_str sprintf('%.4f', row.KDE_corr_NS_BOSS_r) ' & ' sprintf('%.6f', row.KDE_corr_NS_BOSS_p) ' \\\\\n'];
    end
end

latex_str = [latex_str '\\hline\n\\end{tabular}'];

% Save LaTeX table to file
fid = fopen(fullfile(fig_dir, ['circular_statistics_results_' current_participant_name '.tex']), 'w');
fprintf(fid, '%s\n', latex_str);
fclose(fid);
fprintf('LaTeX table saved to %s\n', fullfile(fig_dir, ['circular_statistics_results_' current_participant_name '.tex']));

%% function definitions 
function EEG_lap = prep_ds_dt_lap(EEG, hjorth, down_fs)

if nargin < 2
    hjorth.channel = 'C3';
    hjorth.electrodes = {'FC1', 'CP1', 'FC5', 'CP5'};
end

%% Settings

temp = struct2cell(EEG.chanlocs);
hjorth.channel_idx = find(strcmpi(hjorth.channel, temp(1,:)));
for i = 1:length(hjorth.electrodes)
    hjorth.electrode_idx(i) = find(strcmpi(hjorth.electrodes{i}, temp(1,:)));
end
clear i

%% Downsample

% settings
ds_fs = down_fs;                               % downsampling freqency in Hz

% start timer
tic

% start notification
disp('Start downsampling...')

% Check EEG structure consistency before resampling
EEG = eeg_checkset(EEG);

% downsample data
EEG = pop_resample(EEG, ds_fs);


%% Detrend

% settings
% define prestim period
t.detrend = [min(EEG.times), max(EEG.times)];         % ms time wrt TMS

% start timer
tic

% start notification
disp('Start detrending...')

% detrend data linearly
% To use this: Requires EEGlab with TESA plugin
EEG = tesa_detrend(EEG, 'linear', t.detrend);

% timing notification
fprintf('detrending completed. Expired time: %.0f seconds \n', toc)
%% Laplacian montage

% start timer
tic

% start notification
disp('Start laplacian...')

EEG_lap = EEG;

% Sanity check: EEG.data should be [nbchan x pnts x trials] or [nbchan x pnts]
data = EEG_lap.data;

% If 2D, reshape to 3D for consistent processing
if ndims(data) == 2
    data = reshape(data, size(data,1), size(data,2), 1);
end

% apply laplacian montage
lap_data = data(hjorth.channel_idx,:,:) - ...
    mean(data(hjorth.electrode_idx,:,:), 1);

clear data

% Update EEG_lap structure
EEG_lap.data = lap_data;
EEG_lap.nbchan = 1;
EEG_lap.chanlocs = [hjorth.channel '-centered laplacian'];
EEG_lap.soundLeadfield = [];

% timing notification
fprintf('Laplacian montage completed. Expired time: %.0f seconds \n', toc)


end % eof [EEG_lap] = preprocess(EEG)
