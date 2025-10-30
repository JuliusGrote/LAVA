%% Random session 1 before training 
stimulators.right.arm;
stimulators.right.setAmplitude(SMA_amplitude);
stimulators.left.arm;
stimulators.left.setAmplitude(M1_amplitude);


phase_options=[-2.7489 -1.9635 -1.1781 -0.3927 0.3927 1.1781 1.9635 2.7489];
pair_order=repmat([1:8],1,7);
single_order=repmat([1:8],1,7);
New_pair_order=pair_order(1:50);
New_single_order=single_order(1:50);
single_phase=New_pair_order(randperm(length(pair_order(1:50))));
pair_phase=New_single_order(randperm(length(single_order(1:50))));

Random_index=randperm(100);

All_pusle=[single_phase,pair_phase]
All_phaseID=[ones(1,50),2*ones(1,50)]

New_All_pusle=All_pusle(Random_index);
New_All_phaseID=All_phaseID(Random_index);

savedData_R1=[]
record.emg_R1 = [];
Phase_tolerance = deg2rad(5);
paired_trial_index = 1;
single_trial_index = 1;
savedData_R1.phases = [];
savedData_R1.singlephases = [];
savedData_R1.pairphases = [];
savedData_R1.MEPs.single = [];
savedData_R1.MEPs.paired = [];
savedData_R1.facs = [];
n_trials_paired = 1;
MEPs_P=[]

n_trials_single = 1;%2--0 modified by jiahua 28-07-2022 for testing only pair pusle
% Create stimulus order vector
stim_order = [zeros(1,n_trials_paired) ones(1,n_trials_single)];
stim_order = stim_order(randperm(length(stim_order)));
if(strcmp(bd.armed, 'yes'))
    bd.disarm;
end
for i=1:length(All_phaseID)
    
	    
        if New_All_phaseID(i)==1
            % SINGLE PULSE
			phase_choice = phase_options(All_pusle(i));
            savedData_R1.pairphases = [savedData_R1.pairphases phase_choice];
            notArmed = strcmp(bd.armed, 'no');
            if(notArmed)
                % Re-enable stimulators just in case
                %this.stim_l.arm;
                bd.triggers_remaining = 1;
                bd.configure_time_port_marker([0 1 200]); % TS % 27.04.22 changed the marker from 0 to 200
                bd.alpha.phase_target(1) = phase_choice;
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
                    trial_index = length(record.emg_R1)+1;
                    pause(EMG_pause)
                    try
                        record.emg_R1(trial_index).signal = bd.mep(sc.mep_channel, 100, 100)*5;
                    catch
                        disp('MEP missed. Retrying...')
                        continue
                    end
                    record.emg_R1(trial_index).signal = record.emg_R1(trial_index).signal - mean(record.emg_R1(trial_index).signal(sc.baselineTimeWindow(1) <= sc.emgTimeAxis & sc.emgTimeAxis <= sc.baselineTimeWindow(2)));
                    record.emg_R1(trial_index).type = 1;    % 1 for single pulse
                    temp = extractMEP(record.emg_R1(trial_index).signal,sc.emgTimeAxis,sc.mepTimeWindow,sc.baselineTimeWindow);
                    MEP_S(single_trial_index) = temp.amplitude;
                    single_trial_index = single_trial_index + 1;
                    mep_bug = true;
                    
                    bd.disarm;
                    
                end
            end
            
       
            
            % PAIRED PULSE
        elseif New_All_phaseID(i)==2
            
            phase_choice = phase_options(All_pusle(i));
            savedData_R1.pairphases = [savedData_R1.pairphases phase_choice];
            notArmed = strcmp(bd.armed, 'no');
            if(notArmed)
                %tic
                % Re-enable stimulators just in case
                %stim_l.arm;
                %stim_r.arm;
                bd.triggers_remaining = 1;
                bd.configure_time_port_marker([0 3 2; 0.006 1 1]); % CS TS % 27.04.22 changed the markers from 0,0 to 2,1
                bd.alpha.phase_target(1) = phase_choice;
                bd.alpha.phase_plusminus(1) = Phase_tolerance;
                ITI = sc.iti(1)+ (sc.iti(2)-sc.iti(1)).*rand(1,1);
                tic
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
                    trial_index = length(record.emg_R1)+1;
                    pause(EMG_pause)
                    try
                        record.emg_R1(trial_index).signal = bd.mep(sc.mep_channel, 100, 100)*5;
                    catch
                        disp('MEP missed. Retrying...')
                        continue
                    end
                    record.emg_R1(trial_index).signal = record.emg_R1(trial_index).signal - mean(record.emg_R1(trial_index).signal(sc.baselineTimeWindow(1) <= sc.emgTimeAxis & sc.emgTimeAxis <= sc.baselineTimeWindow(2)));
                    record.emg_R1(trial_index).type = 2;    % 2 for paired pulse
                    temp = extractMEP(record.emg_R1(trial_index).signal,sc.emgTimeAxis,sc.mepTimeWindow,sc.baselineTimeWindow);
                    MEPs_P(paired_trial_index) = temp.amplitude;
                    paired_trial_index = paired_trial_index + 1;
                    mep_bug = true;
                    
                    bd.disarm;
                    
                end
            end
            
              
         end 

%       facilitation_ratio = mean(MEPs_P)/mean(MEP_S);
% %     facilitation_ratio = mean(MEPs_P)/MEP_s(p_id);% single pulse baseline based 
%     savedData.MEPs.single = [savedData.MEPs.single mean(MEP_S)];
%     savedData.MEPs.paired = [savedData.MEPs.paired mean(MEPs_P)];
%     
%     
 disp(i)
end 
   
savedData_R1.All_phaseID=New_All_phaseID;
savedData_R1.All_pusle=New_All_pusle;
savedData_R1.Random_index=Random_index;

savedData_R1.MEPs.single = MEP_S;
savedData_R1.MEPs.paired = MEPs_P;   
savedData_R1.facs = mean(movmean(MEPs_P,3))/mean(movmean(MEP_S,3));

    
    
  
disp('Done')
savedData_Random1=savedData_R1;
saved_dataRandom1 = [folderpath,record.subject,'_',num2str(record.session),'_savedData_Random1','.mat']
save(saved_dataRandom1,'savedData_R1')
