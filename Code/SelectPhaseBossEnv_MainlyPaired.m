classdef SelectPhaseBossEnv_MainlyPaired < rl.env.MATLABEnvironment
    %SELECTPHASEBOSSENV Phase-selection environment wrapping the
    %boss-device; actions are phases (in radians) -- con
    
    properties
        NumberOfSteps;
        InitialPhase = deg2rad(-90);
        RewardCap = 10;
        Phase_tolerance = deg2rad(5);
    end
    
    properties(Access=protected)
        bd;
        stim_r;
        stim_l;
        sc;
        record;
        IsDone = false;
        Trial = 0;
        fac_ratio = [];
        phases = [];
        EMG_pause = 0.4;
        Episode;
        save_flag;
        MEP_S;
        mep_p=[];
        MEPS_all;
        MEPS_all_paired;
        temp_single_MEP;
        Actions;
       reward=[];
       
    end
    
    
    methods
        function this = SelectPhaseBossEnv_MainlyPaired(bd, stim_r, stim_l, sc, record, save_flag, MEP_s,MEPS_all,MEPS_all_paired)
            %SELECTPHASEBOSSENV Wrap given boss-device into a
            %RL-environment
            
           

            
            
            ObservationInfo = rlNumericSpec([2 1]);
            ObservationInfo.Name = 'Observation';
            ObservationInfo.Description = 'phase, facilitation';
            
            % Initialize Action settings
            ActionInfo = rlFiniteSetSpec(1:8);
            ActionInfo.Name = 'Action';
            
            this = this@rl.env.MATLABEnvironment(ObservationInfo,ActionInfo);
            this.bd = bd;
            this.stim_r = stim_r;
            this.stim_l = stim_l;
            this.sc = sc;
            this.record = record;
            this.Episode = 0;
            this.save_flag = save_flag;
            this.MEP_S = MEP_s;
            this.MEPS_all= MEPS_all;
            this.MEPS_all_paired= MEPS_all_paired;
            this.temp_single_MEP =[];
            this.reward=[];
%             this.Action=[];
%             this.Avgfac = Avgfac;
        end
        
        function [Observation,Reward,IsDone,LoggedSignals] = step(this,Action,bd)
            % STEP deliver one paired pulse at phase specified by Action
            % (in radians), and five single pulses
            
            
            
            MEP_P = [];
%             MEP_S = []; commented by jiahua 07-01-2022
            this.Trial = this.Trial + 1;
            LoggedSignals = [];
            phase = (Action - 4.5)*0.25*pi; 
            %phase = deg2rad(Action * 45);
            this.phases = [this.phases phase];
            n_trials_paired = 2;
            n_trials_single = 0;
            % determine if update the single pulse
            if rem(this.Trial,4)==0
                 n_trials_single = 2;
                 
            end
                
            
            % Perform paired pulses
            paired_trial_index = 1;
            single_trial_index = 1;
           
           
            % Create stimulus order vector
            stim_order = [zeros(1,n_trials_paired) ones(1,n_trials_single)];
            stim_order = stim_order(randperm(length(stim_order)));
            
            if(strcmp(this.bd.armed, 'yes'))
                this.bd.disarm;
            end
            
            fprintf('\nStep %d.\n',this.Trial)
            
            temp_single_MEP =[];
            for i = 1:length(stim_order)
                tic
               
                trial_index = [];         
                if stim_order(i)
                    % SINGLE PULSE
                    
                    notArmed = strcmp(this.bd.armed, 'no');
                    if(notArmed)
                        % Re-enable stimulators just in case
                        %this.stim_l.arm;
                        this.bd.triggers_remaining = 1;
                        this.bd.configure_time_port_marker([0 1 200]); % TS % 27.04.22 changed the marker from 0 to 200
                        %this.bd.alpha.phase_target(1) = Action;
                        this.bd.alpha.phase_target(1) = phase; % Note: Dania changed Action here to phase
                        this.bd.alpha.phase_plusminus(1) = this.Phase_tolerance;
                        ITI = this.sc.iti(1)+ (this.sc.iti(2)-this.sc.iti(1)).*rand(1,1);
                        while toc < ITI
                        end
                        %toc % Uncomment to see ITI
                        
                        this.bd.arm;
                        
                    end
                    
                    while(this.bd.triggers_remaining == 1)
                        pause(0.01)
                    end
                    
                    % trigger has been executed, move to the next condition
                    mep_bug = false;
                   
                    while mep_bug == false
                        if(this.bd.triggers_remaining == 0)
                            temp = [];
                            trial_index = length(this.record.emg)+1;
                            pause(this.EMG_pause)
                            try
                                this.record.emg(trial_index).signal = this.bd.mep(this.sc.mep_channel, 100, 100)*5;
                            catch
                                disp('MEP missed. Retrying...')
                                continue
                            end
                            this.record.emg(trial_index).signal = this.record.emg(trial_index).signal - mean(this.record.emg(trial_index).signal(this.sc.baselineTimeWindow(1) <= this.sc.emgTimeAxis & this.sc.emgTimeAxis <= this.sc.baselineTimeWindow(2)));
                            this.record.emg(trial_index).type = 1;    % 1 for single pulse
                            temp = extractMEP(this.record.emg(trial_index).signal,this.sc.emgTimeAxis,this.sc.mepTimeWindow,this.sc.baselineTimeWindow);
                            temp_single_MEP(single_trial_index) = temp.amplitude;
                            single_trial_index = single_trial_index + 1;
                            mep_bug = true;
                            this.bd.disarm;
                            
                        end
                    end
                    
                    % PAIRED PULSE
                else
                    
                    
                    notArmed = strcmp(this.bd.armed, 'no');
                    if(notArmed)
                        %tic
                        % Re-enable stimulators just in case
                        %this.stim_l.arm;
                        %this.stim_r.arm;
                        this.bd.triggers_remaining = 1;
                        this.bd.configure_time_port_marker([0 3 2; 0.006 1 1]); % CS TS % 27.04.22 changed the markers from 0,0 to 2,1
                        %this.bd.alpha.phase_target(1) = Action;
                        this.bd.alpha.phase_target(1) = phase;
                        this.bd.alpha.phase_plusminus(1) = this.Phase_tolerance;
                        ITI = this.sc.iti(1)+ (this.sc.iti(2)-this.sc.iti(1)).*rand(1,1);
                        while toc < ITI
                        end
                        %toc % Uncomment to see ITI
                        
                        this.bd.arm;
                        
                    end
                    
                    while(this.bd.triggers_remaining == 1)
                        pause(0.01)
                    end
                    % trigger has been executed, move to the next
                    % conditionprofile onprofile onprofile on_
                    mep_bug = false;
                    while mep_bug == false
                        if(this.bd.triggers_remaining == 0)
                            temp = [];
                            trial_index = length(this.record.emg)+1;
                            pause(this.EMG_pause)
                            try
                                this.record.emg(trial_index).signal = this.bd.mep(this.sc.mep_channel, 100, 100)*5;
                            catch
                                disp('MEP missed. Retrying...')
                                continue
                            end
                            this.record.emg(trial_index).signal = this.record.emg(trial_index).signal - mean(this.record.emg(trial_index).signal(this.sc.baselineTimeWindow(1) <= this.sc.emgTimeAxis & this.sc.emgTimeAxis <= this.sc.baselineTimeWindow(2)));
                            this.record.emg(trial_index).type = 2;    % 2 for paired pulse
                            temp = extractMEP(this.record.emg(trial_index).signal,this.sc.emgTimeAxis,this.sc.mepTimeWindow,this.sc.baselineTimeWindow);
                            MEP_P(paired_trial_index) = temp.amplitude;
                            paired_trial_index = paired_trial_index + 1;
                            mep_bug = true;
                            
                            this.bd.disarm;
                            
                        end
                    end
                end
            end
             
%             if n_trials_single == 2
%                 this.MEP_S(Action) = mean(movmean([this.MEPS_all temp_single_MEP],3));
%             end
%                 Reward = 0.8*((mean(MEP_P)-1.5*mean(this.MEPS_all_paired))/1000)+0.2*(mean(temp_single_MEP)-mean(this.MEPS_all))/1000;
%                 
%             else 
%                 Reward = (mean(MEP_P)-1.5*mean(this.MEPS_all_paired))/1000+0.2*(this.MEP_S(Action)-mean(this.MEPS_all))/1000;
          
%             
%             fprintf(' Action: %d.\n',  Action)
%             fprintf('temp_single_MEP: %d.\n',temp_single_MEP)
%             fprintf('this.MEP_S(Action): %d.\n', this.MEP_S(Action))
%             fprintf('mean(MEP_P): %d.\n',mean(MEP_P))
            
            
            facilitation_ratio = mean(MEP_P)/this.MEP_S(Action);% modified by jiahua for mainlypaired
            this.mep_p = [this.mep_p mean(MEP_P)];%modified by jiahua for only pair pusle during the traning 
    
            this.fac_ratio = [this.fac_ratio facilitation_ratio];
            this.temp_single_MEP = [this.temp_single_MEP mean(temp_single_MEP)];
            this.Actions = [this.Actions Action];
            
%             Reward = (mean(MEP_P)/1000)*(facilitation_ratio - 1.4725);
%             Reward = mean(MEP_P)/1000;
%             Reward = 0.8.*(mean(MEP_P)/1000)+0.2.*(mean(this.MEPS_all)/1000)-1;
%             Reward = 0.8*(mean(MEP_P)/1000)+0.2*(mean(this.MEP_S)/1000)-1;

%             fprintf(' Mean Pair MEP: %d.\n',  mean(MEP_P))
%             fprintf('Avg pair MEP: %d.\n',mean(this.MEPS_all_paired))
%             fprintf('Single MEP: %d.\n', temp_single_MEP)
%             fprintf('All MEP Avg: %d.\n',mean(this.MEPS_all))
%             Reward =
            Reward =(mean(MEP_P)-1.5*mean(this.MEPS_all_paired))/1000;%enhanced
%           Reward = (0.7*mean(this.MEPS_all_paired)-mean(MEP_P))/1000;%reduced
% % %                 
            
            Reward = sign(Reward) * min(abs(Reward), this.RewardCap);
            this.reward=[this.reward Reward];
            Observation = [phase; facilitation_ratio];
            fprintf('Reward: %d.\n', Reward)
            IsDone = this.Trial >= this.NumberOfSteps;
            this.IsDone = IsDone;
            if(this.IsDone)
                disp('Episode complete.')
                 %%% Initialize boss
            
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
            
            pause(15)
            bd.arm;
            
            end
        end
        
        % Reset environment to initial state and output initial observation
        function InitialObservation = reset(this)
            this.Trial = 0;
            this.IsDone = false;
            this.Episode = this.Episode + 1;
            InitialObservation = [this.InitialPhase; 1];
        end
        
        function saveData = saveFacRatios_Phases(this)
            this.record.fac_rat = this.fac_ratio;
            this.record.phases = this.phases;
            this.record.temp_single_MEP = this.temp_single_MEP;
            this.record.Actions = this.Actions;
            this.record.reward=this.reward;
            % Save emg data
            if this.save_flag
                if ~exist([pwd '\saved'],'dir')
                    mkdir('saved')
                end
                filename = fullfile(pwd, '\saved\',['CHASMEC_', num2str(this.record.subjectID),'\','Session_' ,num2str(this.record.session)],'\',[this.record.subject '_session' num2str(this.record.session) '.mat']);
                saveData = this.record;
                save(filename,'saveData')
                disp('Record saved.')
            end
        end
        
        function setNumStepsPerEp(this, numSteps)
            this.NumberOfSteps = numSteps;
        end
    end
end