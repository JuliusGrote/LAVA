
%% Initialize boss
bd = bossdevice;
bd.eeg_channels = 64;
bd.aux_channels = 1;
%set up spatial filter to detect the mu-rythm
spf = zeros(64, 1);
% [C3 FC1 CP1 FC5 CP5] (C3 Hjorth filter)
spf([5 21 23 25 27], 1) = [1 -0.25 -0.25 -0.25 -0.25]; % left
% [C4 FC2 CP2 FC6 CP6] (C4 Hjorth filter)
%spf([6 22 24 26 28], 2) = [1 -0.25 -0.25 -0.25 -0.25]; % right
bd.spatial_filter_weights = spf;

bd.alpha.offset_samples = 6; %depends on loop delay, they have demo script for this, I keep 6 from Gabor
bd.alpha.amplitude_min(1) = 0;
bd.alpha.amplitude_max(1) = 1e6;
bd.theta.ignore
bd.beta.ignore
bd.armed = 'no';
bd.min_inter_trig_interval = 1.5;
phase_tolerance = pi/40;
bd.alpha.phase_plusminus(1) = phase_tolerance;
%bd.sample_and_hold_period = 0.1;

sc=[]; %pool all properties we cant put into bd here
% sc.iti = [0.75 1.25]; % [4 6]: ITI limits (s)  % inter trial interval
sc.iti = [1 2]; % [4 6]: ITI limits (s)  % inter trial interval
sc.samplingFrequency = 1000; % 5000 sampling frequency (Hz)
sc.emgTimeWindow = [-100 100]; % [-100 100] time window of the EMG trial (ms)
sc.mepTimeWindow = [18 55]; % [18 55] time window for MEP detection (ms)
sc.baselineTimeWindow = [-100 0]; % [-100 0] time window for baseline analysis (ms)
sc.baselineThreshold = 15; % 15 threshold absolute EMG amplitude (?V) in the baseline time window
sc.applyBaselineThreshold = 1; % 1 remeasure trials that have bad baseline (0 = do not apply this threshold)
sc.emgTimeAxis = createTimeAxis(sc.emgTimeWindow,sc.samplingFrequency); % time axis for EMG data (ms)

sc.mep_channel = 1; %<--- SET
%% Random session3 after traning 
stimulators.right.arm;
stimulators.right.setAmplitude(SMA_amplitude);
stimulators.left.arm;
stimulators.left.setAmplitude(M1_amplitude);


phase_options=[-2.7489 -1.9635 -1.1781 -0.3927 0.3927 1.1781 1.9635 2.7489];


pair_order=repmat([1:8],1,50);

pair_phase=pair_order(randperm(length(pair_order)));


savedData_R_training=[];
record.emg_R_training = [];
Phase_tolerance = deg2rad(5);
paired_trial_index = 1;
single_trial_index = 1;

savedData_R_training.singlephases = [];
savedData_R_training.pairphases = [];
savedData_R_training.MEPs.single = [];
savedData_R_training.MEPs.paired = [];


MEPs_PR=[];
MEPs_SR=[];

pairphases=[];
singlephases =[];





% Perform paired pulses
paired_trial_index = 1;
single_trial_index = 1;


if(strcmp(bd.armed, 'yes'))
    bd.disarm;
end
%% run the main random part 
for j=220:length(pair_phase)
% for j=1:5
    
      %determine if update the single pulse
            if rem(j,4)==0
                 n_trials_single = 2;
                 n_trials_paired = 2;
            else
               n_trials_single = 0;
               n_trials_paired = 2;
            end
                
%             n_trials_single = 0;
%             n_trials_paired = 2;
           
            phasetarget=phase_options(pair_phase(j));
           
            % Create stimulus order vector
            stim_order = [zeros(1,n_trials_paired) ones(1,n_trials_single)];
            stim_order = stim_order(randperm(length(stim_order)));
            
            if(strcmp(bd.armed, 'yes'))
                bd.disarm;
            end
            
            fprintf('\nStep %d.\n',j)
        
            for i = 1:length(stim_order)
                tic
               
                trial_index = [];         
                if stim_order(i)
                    % SINGLE PULSE
                    
                    notArmed = strcmp(bd.armed, 'no');
                    if(notArmed)
                        % Re-enable stimulators just in case
                        %this.stim_l.arm;
                        bd.triggers_remaining = 1;
                        bd.configure_time_port_marker([0 1 200]); % TS % 27.04.22 changed the marker from 0 to 200
                        %this.bd.alpha.phase_target(1) = Action;
                        bd.alpha.phase_target(1) = phasetarget; % Note: Dania changed Action here to phase
                        bd.alpha.phase_plusminus(1) = Phase_tolerance;
                        ITI = sc.iti(1)+ (sc.iti(2)-sc.iti(1)).*rand(1,1);
                        while toc < ITI
                        end
                        %toc % Uncomment to see ITI
                        
                        bd.arm;
                        
                    end
                    
                    while(bd.triggers_remaining == 1)
                        pause(0.01)
                    end
                    
                    % trigger has been executed, move to the next condition
                    mep_bug = false;
                   
                    while mep_bug == false
                        if(bd.triggers_remaining == 0)
                            temp = [];
                            trial_index = length(record.emg_R_training)+1;
                            pause(EMG_pause)
                            try
                               record.emg_R_training(trial_index).signal = bd.mep(sc.mep_channel, 100, 100)*5;
                            catch
                                disp('MEP missed. Retrying...')
                                continue
                            end
                            record.emg_R_training(trial_index).signal = record.emg_R_training(trial_index).signal - mean(record.emg_R_training(trial_index).signal(sc.baselineTimeWindow(1) <= sc.emgTimeAxis & sc.emgTimeAxis <= sc.baselineTimeWindow(2)));
                            record.emg_R_training(trial_index).type = 1;    % 1 for single pulse
                            temp = extractMEP(record.emg_R_training(trial_index).signal,sc.emgTimeAxis,sc.mepTimeWindow,sc.baselineTimeWindow);
                            MEPs_SR(single_trial_index) = temp.amplitude;
                            singlephases=[singlephases,phasetarget];
                            single_trial_index = single_trial_index + 1;
                            mep_bug = true;
                            bd.disarm;
                            
                        end
                    end
                    
                    % PAIRED PULSE
                else
                    
                    
                    notArmed = strcmp(bd.armed, 'no');
                    if(notArmed)
                        %tic
                        % Re-enable stimulators just in case
                        %this.stim_l.arm;
                        %this.stim_r.arm;
                        bd.triggers_remaining = 1;
                        bd.configure_time_port_marker([0 3 2; 0.006 1 1]); % CS TS % 27.04.22 changed the markers from 0,0 to 2,1
                        %this.bd.alpha.phase_target(1) = Action;
                        bd.alpha.phase_target(1) = phasetarget;
                        bd.alpha.phase_plusminus(1) =Phase_tolerance;
                        ITI = sc.iti(1)+ (sc.iti(2)-sc.iti(1)).*rand(1,1);
                        while toc < ITI
                        end
                        %toc % Uncomment to see ITI
                        
                      bd.arm;
                        
                    end
                    
                    while(bd.triggers_remaining == 1)
                        pause(0.01)
                    end
                    % trigger has been executed, move to the next condition
                    mep_bug = false;
                    while mep_bug == false
                        if(bd.triggers_remaining == 0)
                            temp = [];
                            trial_index = length(record.emg_R_training)+1;
                            pause(EMG_pause)
                            try
                               record.emg_R_training(trial_index).signal = bd.mep(sc.mep_channel, 100, 100)*5;
                            catch
                                disp('MEP missed. Retrying...')
                                continue
                            end
                            record.emg_R_training(trial_index).signal = record.emg_R_training(trial_index).signal - mean(record.emg_R_training(trial_index).signal(sc.baselineTimeWindow(1) <= sc.emgTimeAxis & sc.emgTimeAxis <= sc.baselineTimeWindow(2)));
                            record.emg_R_training(trial_index).type = 2;    % 2 for paired pulse
                            temp = extractMEP(record.emg_R_training(trial_index).signal,sc.emgTimeAxis,sc.mepTimeWindow,sc.baselineTimeWindow);
                            MEPs_PR(paired_trial_index) = temp.amplitude;
                            pairphases=[pairphases,phasetarget];
                            paired_trial_index = paired_trial_index + 1;
                            mep_bug = true;
                            
                            bd.disarm;
                            
                        end
                    end
                end
            end
   
end 
 

% save all the data 
savedData_R_training.singlephases = singlephases;
savedData_R_training.pairphases = pairphases;
savedData_R_training.MEPs.single = MEPs_SR;
savedData_R_training.MEPs.paired = MEPs_PR;

savedData_Random_training=savedData_R_training;
savedData_Random_training = [folderpath,record.subject,'_',num2str(record.session),'_savedData_Random_traning','.mat'];
save(savedData_Random_training,'savedData_R_training')
