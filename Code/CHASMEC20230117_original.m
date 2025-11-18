% 2022-10-26 CHASMEC with bossdevice v2021a

 cd C:\Users\bnp-admin\Desktop\Projects\CHASMEC
 delete(instrfindall) % use this commmand in begiinning


if (~exist('stimulators') | isempty('stimulators'))
    try
        stimulators = [];
        stimulators.left = magventure('COM3');
        stimulators.right = magventure('COM4');
        
        stimulators.left.connect();
        stimulators.right.connect();
%         stimulators.right.setWaveform('Monophasic','Normal');
    
    catch
        fprintf(2, 'Problem with the MAGIC toolbox, please restart Matlab.\n');
    end
end

%% must check the pulse type.
% stimulators.right.setWaveform('Monophasic','00')
% stimulators.right.setWaveform('Monophasic','Normal',1, 1, 0.2,1)
% stimulators.right.setPage('Main')

stimulators.right.setPage('Main')
stimulators.right.setWaveform('Monophasic','Normal',2,2,1,0)


% stimulators.right.getWaveform()
% must be monophasic for SMA 

%% subject ID 
rng(now)

save_flag = 1;

record = [];

% ---------- below set by user ----------
subject_number =47; %<--- SET
session_number = 2; %<--- SET

%%  save the reward function under this session 
% ---------- above set by user ----------
%           R1  Reward = facilitation_ratio - 1.5;
%           R2  Reward = (mean(MEP_P)/1000)*(facilitation_ratio - 1.4725);
%           R3  Reward = mean(MEP_P)/1000;
%           R4  Reward = 0.8.*(mean(MEP_P)/1000)+0.2.*(mean(this.MEPS_all)/1000)-1;
%           R5  Reward = 0.8*(mean(MEP_P)/1000)+0.2*(mean(this.MEP_S)/1000)-1;
%           R6  Reward = 0.8*((mean(MEP_P)-1.5*mean(this.MEPS_all_paired))/1000)+0.2*(mean(temp_single_MEP)-mean(this.MEPS_all))/1000;
%           R7  Reward = (mean(MEP_P)-1.5*mean(this.save(fullfile(subject_folder, 'trigger_times.mat'), 'trigger_timestamps');MEPS_all_paired))/1000
%           R8  Reward = (0.7*mean(this.MEPS_all_paired)-(mean(MEP_P))/1000
%           R9  Random

record.condition = 'R7'; % 
record.IsExperiment = 'Large experiment trials'; % 
record.pusle='2 pair per step , only 2 single after 4steps';
record.emg = [];
record.seq='Enhanced_';% reduce, random , we need 2 enhance and 1 reduce and 1 random session 
%% make folder for all matlab data under this subject 
if ~exist([pwd '\saved\','CHASMEC_' ,num2str(subject_number),'\Session_',num2str(session_number)],'dir')
    mkdir([pwd '\saved\','CHASMEC_' ,num2str(subject_number),'\Session_',num2str(session_number),'\'])
end

% make folder in the data server for neurone data 

if ~exist(['Z:\Experimental Data\2021-10 CHASMEC\NeuroOne data\','CHASMEC_' ,num2str(subject_number),'\',record.seq ,num2str(session_number)],'dir')
    mkdir(['Z:\Experimental Data\2021-10 CHASMEC\NeuroOne data\','CHASMEC_' ,num2str(subject_number),'\',record.seq ,num2str(session_number),'\Resting\'])
    mkdir(['Z:\Experimental Data\2021-10 CHASMEC\NeuroOne data\','CHASMEC_' ,num2str(subject_number),'\',record.seq ,num2str(session_number),'\Training\'])
    mkdir(['Z:\Experimental Data\2021-10 CHASMEC\NeuroOne data\','CHASMEC_' ,num2str(subject_number),'\',record.seq ,num2str(session_number),'\PI1\'])
    mkdir(['Z:\Experimental Data\2021-10 CHASMEC\NeuroOne data\','CHASMEC_' ,num2str(subject_number),'\',record.seq ,num2str(session_number),'\PI2\'])
end


folderpath=[pwd '\saved\','CHASMEC_' ,num2str(subject_number),'\Session_',num2str(session_number),'\']; 

record.subject = ['CHASMEC_' , num2str(subject_number,'%03.f')];
record.session = session_number;
record.subjectID=subject_number;



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
sc.iti = [1 2]; % [4 6]: ITI limits (s)  % inter trial interval, 1-2 setting demos 2-3 seconds in practice as jetter, code runing delay, boss delay make it round 2-3 seconds 
sc.samplingFrequency = 1000; % 5000 sampling frequency (Hz)
sc.emgTimeWindow = [-100 100]; % [-100 100] time window of the EMG trial (ms)
sc.mepTimeWindow = [18 55]; % [18 55] time window for MEP detection (ms)
sc.baselineTimeWindow = [-100 0]; % [-100 0] time window for baseline analysis (ms)
sc.baselineThreshold = 15; % 15 threshold absolute EMG amplitude (?V) in the baseline time window
sc.applyBaselineThreshold = 1; % 1 remeasure trials that have bad baseline (0 = do not apply this threshold)
sc.emgTimeAxis = createTimeAxis(sc.emgTimeWindow,sc.samplingFrequency); % time axis for EMG data (ms)

sc.mep_channel = 1; %<--- SET

%% Add variables
EMG_pause = 0.3;   % Pause after trigger to save EMG epoch
fac_ratio = [];
n_single_blocks = 1;
n_paired_blocks = 1;
MEPs = [];
sc.single_data.emg=[];
sc.paired_data.emg=[];
MEP_s=[];

%% Find motor hotspot
% Note: NeurOne needs to be running as the triggers pass through the EEG system
bd.configure_time_port_marker([0 1 0]); % TS only NOTE: Changed from 1 to 2 by Dania, 14.03.22
stimulators.left.arm;
stimulators.left.setAmplitude(65);


while(true), bd.manualTrigger, fprintf('.\n'), pause(2), end

% 
% % Note: NeurOne n eeds to be running as the triggers pass through the EEG system
% bd.configure_time_port_marker([0 3 0]); % TS only NOTE: Changed from 1 to 2 by Dania, 14.03.22
% stimulators.right.arm;
% stimulators.right.setAmplitude(75);
% 
% 
% while(true), bd.manualTrigger, fprintf('.\n'), pause(2), end

%% Find SMA hotspot
M1_amplitude=70; % SET
SMA_amplitude=70; % SET

% Ready stimulators
bd.configure_time_port_marker([0 3 2; 0.006 1 1]); % CS TS % 27.04.22 changed the markers from 0,0 to 2,1
stimulators.right.arm;
stimulators.right.setAmplitude(SMA_amplitude);
stimulators.left.arm;
stimulators.left.setAmplitude(M1_amplitude);
    
while(true), bd.manualTrigger, fprintf('.\n'), pause(2), end

%% Clear stuff
clear sc.single_data

%% M1 motor threshold

M1_amplitude=75;    % SET
trials_single = 10;
MT = 1000; % micro volts

% Ready stimulator 
bd.configure_time_port_marker([0 1 200]); % TS only NOTE: Changed from 1 to 2 by Dania, 14.03.22

% To test the phase after the intervention,comment this after use:
% this.bd.alpha.phase_target(1) = 2.7489;
% this.bd.alpha.phase_plusminus(1) = deg2rad(5);

for s_tr=1:trials_single
    stimulators.left.arm;
    stimulators.left.setAmplitude(M1_amplitude);
    
    bd.manualTrigger
    pause(EMG_pause);
    sc.single_data(n_single_blocks).emg(s_tr,:) = bd.mep(sc.mep_channel, 100, 100)*5; % Read MEP
    sc.single_data(n_single_blocks).emg(s_tr,:) = sc.single_data(n_single_blocks).emg(s_tr,:) ...
        - mean(sc.single_data(n_single_blocks).emg(s_tr,sc.baselineTimeWindow(1) <= sc.emgTimeAxis & sc.emgTimeAxis <= sc.baselineTimeWindow(2)));
    sc.single_data(n_single_blocks).mep(s_tr) = extractMEP(sc.single_data(n_single_blocks).emg(s_tr,:),sc.emgTimeAxis,sc.mepTimeWindow,sc.baselineTimeWindow);
    
    rand_ITI = sc.iti(1)+ (sc.iti(2)-sc.iti(1)).*rand(1,1); %Assigning New Random ITI for this Trial to the BOSS Device
    pause(rand_ITI)
end

% MEP size data
MEPs.single(n_single_blocks).samples = [sc.single_data(n_single_blocks).mep.amplitude];
MEPs.single(n_single_blocks).avg = mean([sc.single_data(n_single_blocks).mep.amplitude]);
MEPs.single(n_single_blocks).med = median([sc.single_data(n_single_blocks).mep.amplitude]);

%Plot
figure(1); clf; hold on
plot(1:n_single_blocks,[MEPs.single.avg],'or')
plot(1:n_single_blocks,[MEPs.single.med],'xb')
plot(1:n_single_blocks,reshape([MEPs.single.samples],trials_single,n_single_blocks),'.black')
xlim([0 n_single_blocks + 1]);
title('MEP amplitudes with single pulse')
ylabel('MEP amplitude (\muV)')
xlabel('Single pulse trial')
legend('Average MEP size','Median MEP size','Single responses')

%plot(sc.single_data(n_single_blocks).emg)
n_over_MT = length(MEPs.single(n_single_blocks).samples(MEPs.single(n_single_blocks).samples > MT));

n_single_blocks = n_single_blocks +1;

fprintf(['%d/%d over %d ' char(181),'V\n'],n_over_MT,trials_single,MT)

saveas(gcf,[folderpath,record.subject,num2str(record.session),'_RMT and AMT for single pusle.png'])

%% SMA-M1 paired pulse
% The SMA+M1 mysteriously saves double MEPs...
% Can it be that mep() uses both markers to get MEPs?

M1_amplitude=75; % SET
SMA_amplitude=81; % SET
trials_paired = 10;

% Ready stimulators
bd.configure_time_port_marker([0 3 2; 0.006 1 1]); % CS TS % 27.04.22 changed the markers from 0,0 to 2,1
    stimulators.right.arm;
    stimulators.right.setAmplitude(SMA_amplitude);
    stimulators.left.arm;
    stimulators.left.setAmplitude(M1_amplitude);

% To test the phase after the intervention:
% this.bd.alpha.phase_target(1) = 2.7489;
% this.bd.alpha.phase_plusminus(1) = deg2rad(5);



tic;
for p_tr=1:trials_paired
    isiList(p_tr) = rand*(sc.iti(2)-sc.iti(1))+sc.iti(1);
    while toc < isiList(p_tr)
    end
    % administer a pulse
    bd.manualTrigger;
    pause(EMG_pause)
    % reset ISI counter
    tic;
    sc.paired_data(n_paired_blocks).emg(p_tr,:) = bd.mep(sc.mep_channel, 100, 100)*5;
    sc.paired_data(n_paired_blocks).emg(p_tr,:) = sc.paired_data(n_paired_blocks).emg(p_tr,:) - mean(sc.paired_data(n_paired_blocks).emg(p_tr,sc.baselineTimeWindow(1) <= sc.emgTimeAxis & sc.emgTimeAxis <= sc.baselineTimeWindow(2)));
    sc.paired_data(n_paired_blocks).mep(p_tr) = extractMEP(sc.paired_data(n_paired_blocks).emg(p_tr,:),sc.emgTimeAxis,sc.mepTimeWindow,sc.baselineTimeWindow);
end

% MEP size data
MEPs.paired(n_paired_blocks).samples=[sc.paired_data(n_paired_blocks).mep.amplitude];
MEPs.paired(n_paired_blocks).avg = mean([sc.paired_data(n_paired_blocks).mep.amplitude]);
MEPs.paired(n_paired_blocks).med = median([sc.paired_data(n_paired_blocks).mep.amplitude]);

% SMA position search with triggered pulses, observe EMG
% Abort using Ctrl+C
%fprintf('\nHotspot search. Press Ctrl+C to abort\n')
%while(true),  bd.manualTrigger(), fprintf('.'), pause(3), end

% Plot
figure(2);clf;hold on
plot(1:n_paired_blocks,[MEPs.paired.avg],'or')
plot(1:n_paired_blocks,[MEPs.paired.med],'xb')
plot(1:n_paired_blocks,reshape([MEPs.paired.samples],trials_paired,n_paired_blocks),'.black')
xlim([0 n_paired_blocks + 1]);
title('MEP amplitudes with paired pulse')
ylabel('MEP amplitude (\muV)')
xlabel('Paired pulse trial')
legend('Average MEP size','Median MEP size','Single responses')

saveas(gcf,[folderpath,record.subject,num2str(record.session),'_pairpusle checking before traning for single pusle.png'])


figure(3);clf;hold on
if(n_single_blocks > 1)
    fac_ratio(n_paired_blocks) = MEPs.paired(n_paired_blocks).avg/MEPs.single(n_single_blocks-1).avg;
    plot(fac_ratio,'xb')
    title('Facilitation ratio')
    ylabel('Relative mean MEP size')
    xlabel('Paired pulse trial')
    xlim([0 n_paired_blocks + 1])
else
    fac_ratio(n_paired_blocks) = 0;
end

n_paired_blocks = n_paired_blocks + 1;

saveas(gcf,[folderpath,record.subject,num2str(record.session),'_fac_checking before traning for single pusle.png'])


%% Random session 1 before training 
Random1
%% Random session 2 before training  
Random2

%% visualize the random data 
single_index1=savedData_R1.All_pusle(savedData_R1.All_phaseID==1)
pair_index1=savedData_R1.All_pusle(savedData_R1.All_phaseID==2)

single_index2=savedData_R2.All_pusle(savedData_R2.All_phaseID==1)
pair_index2=savedData_R2.All_pusle(savedData_R2.All_phaseID==2)

figure(56)

SPallr{1,1}=savedData_R1.MEPs.paired';
SPallr{1,2}=savedData_R2.MEPs.paired';



SPallr{2,1}=savedData_R1.MEPs.single';
SPallr{2,2}=+savedData_R2.MEPs.single';



xlab={'R1','R2'};
col=[
102,194,165,200;
141,160,203,200];
col=col/255;

multiple_boxplot(SPallr',xlab,{'Single pulse','Pair pulse'},col')
xlabel('Phase')
ylabel('MEPs')

saveas(gcf,[folderpath, num2str(record.subject),num2str(record.session),'_Single VS paried MEPS AVG for the baseline only.png'])




%% Calulate the avg single mep per phase and prepare all the single and pair MEP to reward function in case 


MEP_s=[]

for j=1:8
    
       MEP_Single_index1=savedData_Random1.All_pusle(savedData_Random1.All_phaseID==1)  
       MEP_Single_index2=savedData_Random2.All_pusle(savedData_Random2.All_phaseID==1)
       MEP_s1(j)=median(movmedian(savedData_Random1.MEPs.single(MEP_Single_index1==j),3))
       MEP_s2(j)=median(movmedian(savedData_Random2.MEPs.single(MEP_Single_index2==j),3))
       
end 
 
 % Avg single pusle 
 MEP_s=(MEP_s1+MEP_s2)/2
 % All the mep of single and paired for reward function 
 MEPS_all=[savedData_Random1.MEPs.single savedData_Random2.MEPs.single];
 MEPS_all_paired=[savedData_Random1.MEPs.paired savedData_Random2.MEPs.paired];
 
 
 
 

MEP_PP=[]

for j=1:8
    
       MEP_p_index1=savedData_R1.All_pusle(savedData_R1.All_phaseID==2)  
       MEP_p_index2=savedData_R2.All_pusle(savedData_R2.All_phaseID==2)
       MEP_sp1(j)=median(movmedian(savedData_R1.MEPs.paired(MEP_p_index1==j),3))
       MEP_sp2(j)=median(movmedian(savedData_R2.MEPs.paired(MEP_p_index2==j),3))
       
end 
 
 % Avg pair pusle 
 MEP_PP=(MEP_sp1+MEP_sp2)/2
 

%% RL Environment preparation

oriMEP=MEP_s;

addpath('BOSS_ML');

%env = SelectPhaseBossEnvBugCatch(bd,stimulators.right,stimulators.left,sc,record,save_flag);
%env = SelectPhaseBossEnv(bd,stimulators.right,stimulators.left,sc,record,save_flag);
env = SelectPhaseBossEnv_MainlyPaired(bd,stimulators.right,stimulators.left,sc,record,save_flag,MEP_s,MEPS_all,MEPS_all_paired);
% env = SelectPhaseBossEnv_Paironly(bd,stimulators.right,stimulators.left,sc,record,save_flag,MEP_s);

actionDef = getActionInfo(env);
observationDef = getObservationInfo(env);

%episodeLength = 200;
%env.EpisodeLength = episodeLength;
%% RL agent preparation
criticNW = [
    sequenceInputLayer(observationDef.Dimension(1), 'Normalization', 'none', 'Name', 'state')
    fullyConnectedLayer(12, 'Name', 'CriticStateFC1')
    leakyReluLayer(0.3, 'Name','Critic Leaky Relu1')
    fullyConnectedLayer(8, 'Name', 'CriticStateFC2')
    leakyReluLayer(0.3, 'Name','Critic Leaky Relu2')
    lstmLayer(10,'OutputMode','sequence','Name','CriticLSTM'); % previously: 10 units
    fullyConnectedLayer(8,'Name','Critic FC3')
    leakyReluLayer(0.3, 'Name','Critic Leaky Relu3')
    fullyConnectedLayer(numel(actionDef.Elements),'Name','output')
    ];

criticOptions = rlRepresentationOptions('LearnRate', 1e-3, ...
    'GradientThreshold', 1, 'UseDevice', 'cpu');
critic = rlQValueRepresentation(criticNW, observationDef, actionDef, ...
                                'Observation', 'state', criticOptions);

agentOptions = rlDQNAgentOptions(...
    'UseDoubleDQN',false, ...
    'TargetSmoothFactor',5e-3, ...
    'ExperienceBufferLength',50, ...
    'MiniBatchSize',32 , ...
    'ResetExperienceBufferBeforeTraining', false, ...
    'SequenceLength',2); % <- this is required!
agentOptions.EpsilonGreedyExploration.EpsilonDecay = 0.005;
rlDQNAgentOptions('DiscountFactor',1.0)
agent = rlDQNAgent(critic, agentOptions);

%% Agent training options

trainOpts = rlTrainingOptions;

trainOpts.MaxEpisodes = 50;
trainOpts.MaxStepsPerEpisode = 40;
trainOpts.StopTrainingCriteria = "EpisodeCount";
trainOpts.StopTrainingValue =10; %10;
trainOpts.ScoreAveragingWindowLength = 1;
trainOpts.Verbose = false;
trainOpts.Plots = "training-progress";
trainOpts.ParallelizationOptions.Mode = "async";

env.setNumStepsPerEp(trainOpts.MaxStepsPerEpisode);


%% Start experiment
% profile on
% Ready stimulators
stimulators.right.arm;
stimulators.right.setAmplitude(SMA_amplitude);
stimulators.left.arm;
stimulators.left.setAmplitude(M1_amplitude);

% Loading the agent and starting the training %TODO: include both options
% to load or create new

agent_name = [folderpath,'dqn_agent_', record.subject,num2str(record.session),'.mat']

if exist(agent_name,'file')
    fprintf('loading agent')
    load(agent_name,'agent')    
end
trainStats = train(agent, env, trainOpts);

if save_flag
    savedData = env.saveFacRatios_Phases();
   
    savedData.record = record; % fac rand (inh)
    
    save(agent_name,'agent')
    Save_name_data = [folderpath,record.subject,'_',num2str(record.session),'_savedData','.mat']
    save(Save_name_data,'savedData')  
    
end

trainStats.agentOptions=agentOptions;
trainStats.MEP_s=MEP_s;
trainStats.MEPsbeforetraning=oriMEP;
trainStats.fac_ratio=fac_ratio;
Save_name = [folderpath,record.subject,'_',num2str(record.session),'trainStats.mat']%% in case overwrite the saved mat file
if exist(Save_name,'file')
  Save_name = [folderpath,record.subject,'_',num2str(record.session),'_',num2str(rand(1)),'trainStats.mat']
  save(Save_name,'trainStats')  
else 
  save(Save_name,'trainStats')  
end 


%% Plot the histograms for the results

% figure(4);
% polarhistogram(savedData.phases,'BinEdges',[-2.7489 -1.9635 -1.1781 -0.3927 0.3927 1.1781 1.9635 2.7489])
% title('Number of pulses per phase')
% saveas(gcf,[pwd '\saved\','dqn_agent_', record.subject,num2str(record.session),'_polarhist.png'])
% 
% figure(5);
% histogram(rad2deg(savedData.phases),'BinEdges',rad2deg([-2.7489 -1.9635 -1.1781 -0.3927 0.3927 1.1781 1.9635 2.7489]))
% title('Number of pulses per phase')
% saveas(gcf,[pwd '\saved\','dqn_agent_', record.subject,num2str(record.session),'_Histogram.png'])
% 
% figure(5);
% histogram(rad2deg(savedData.phases))
% title('Number of pulses per phase')
% saveas(gcf,[pwd '\saved\','dqn_agent_', record.subject,num2str(record.session),'_Histogram.png'])
% 
% 
% 
% figure(6);
% histogram(savedData.phases(98:121),'BinEdges',[-2.7489 -1.9635 -1.1781 -0.3927 0.3927 1.1781 1.9635 2.7489])

%% plot total nums of each phases
EP=10;
stepnum=40;

figure(7);
[a,b]=hist(savedData.phases,unique(savedData.phases))


X = categorical({'-157.5^\circ', '-112.5^\circ',  '-67.5^\circ',  '-22.5^\circ',   '22.5^\circ',   '67.5^\circ',  '112.5^\circ',  '157.5^\circ'});
X = reordercats(X,{'-157.5^\circ', '-112.5^\circ',  '-67.5^\circ',  '-22.5^\circ',   '22.5^\circ',   '67.5^\circ',  '112.5^\circ',  '157.5^\circ'});
bar(X,a')
% legend({'EP1','EP2','EP3','EP4','EP5'})
xlabel('Phase');
ylabel('Nums');
title('Number of pulses per phase')
saveas(gcf,[folderpath, record.subject,num2str(record.session),'_Histogram_total_nums_of_each_phase23333.png'])



%% check the optimized value per episode 
 
for EPid=1:10
[a,b]=hist(savedData.phases((EPid-1)*40+1:EPid*40),unique(savedData.phases))


X = categorical({'-157.5^\circ', '-112.5^\circ',  '-67.5^\circ',  '-22.5^\circ',   '22.5^\circ',   '67.5^\circ',  '112.5^\circ',  '157.5^\circ'});
X = reordercats(X,{'-157.5^\circ', '-112.5^\circ',  '-67.5^\circ',  '-22.5^\circ',   '22.5^\circ',   '67.5^\circ',  '112.5^\circ',  '157.5^\circ'});
% bar(X,a')
bar(a)
% legend({'EP1','EP2','EP3','EP4','EP5'})
xlabel('Phase');
ylabel('Nums');
title(['Number of pulses per phase in EP\_',num2str(EPid)])
saveas(gcf,[folderpath,num2str(record.subject),'_',num2str(record.session),'EP',num2str(EPid),'_Num_per_Phase_Histogram1.png'])
close all 
end 


%% % plot the fac of the phases in EP
figure(8);
fac=savedData.fac_rat;
phase=rad2deg(savedData.phases);

epallm=[];
phaselabel={'-157.5^\circ', '-112.5^\circ',  '-67.5^\circ',  '-22.5^\circ',   '22.5^\circ',   '67.5^\circ',  '112.5^\circ',  '157.5^\circ'};



for i=1:EP
    Fstart=stepnum*(i-1)+1;
    Fend=stepnum*i;
    Fac_sub=fac(Fstart:Fend);
    epall{i,1}=Fac_sub(find(phase(Fstart:Fend)==-157.5000));
    epall{i,2}=Fac_sub(find(phase(Fstart:Fend)==-112.5000));
    epall{i,3}=Fac_sub(find(phase(Fstart:Fend)==-67.5000));
    epall{i,4}=Fac_sub(find(phase(Fstart:Fend)==-22.5000));
    epall{i,5}=Fac_sub(find(phase(Fstart:Fend)==22.5000));
    epall{i,6}=Fac_sub(find(phase(Fstart:Fend)==67.5000));
    epall{i,7}=Fac_sub(find(phase(Fstart:Fend)==112.5000));
    epall{i,8}=Fac_sub(find(phase(Fstart:Fend)==157.5000));


    epallm(i,1)=median(Fac_sub(find(phase(Fstart:Fend)==-157.5000)));
    epallm(i,2)=median(Fac_sub(find(phase(Fstart:Fend)==-112.5000)));
    epallm(i,3)=median(Fac_sub(find(phase(Fstart:Fend)==-67.5000)));
    epallm(i,4)=median(Fac_sub(find(phase(Fstart:Fend)==-22.5000)));
    epallm(i,5)=median(Fac_sub(find(phase(Fstart:Fend)==22.5000)));
    epallm(i,6)=median(Fac_sub(find(phase(Fstart:Fend)==67.5000)));
    epallm(i,7)=median(Fac_sub(find(phase(Fstart:Fend)==112.5000)));
    epallm(i,8)=median(Fac_sub(find(phase(Fstart:Fend)==157.5000)));

end 


X = categorical({'-157.5^\circ', '-112.5^\circ',  '-67.5^\circ',  '-22.5^\circ',   '22.5^\circ',   '67.5^\circ',  '112.5^\circ',  '157.5^\circ'});
X = reordercats(X,{'-157.5^\circ', '-112.5^\circ',  '-67.5^\circ',  '-22.5^\circ',   '22.5^\circ',   '67.5^\circ',  '112.5^\circ',  '157.5^\circ'});
bar(X,epallm')
legend(['Episode 1-',num2str(EP)])
xlabel('Phase');
ylabel('Facilitation');

title('Fac of pulses per phase')
saveas(gcf,[folderpath, record.subject,num2str(record.session),'_Histogram_fac_of_each_phase.png'])

%% plot the fac of the phases in EP in boxplot 
figure(9)
phases = [-157.5, -112.5, -67.5, -22.5, 22.5, 67.5, 112.5, 157.5];
phase_label = [];
num = [];
data = [];
% EP=4;
% stepnum=40;
for i = 1:EP
    Fstart=stepnum*(i-1)+1;
    Fend=stepnum*i;
    Fac_sub=fac(Fstart:Fend);
    for j = 1:8
        temp = Fac_sub(find(phase(Fstart:Fend)==phases(j)));
        % Create a phase label for each data point found
        phase_label = [phase_label,ones(1,length(temp))*phases(j)];
        % create an index label for each data point found
        num = [num,ones(1,length(temp))*i];
        % save the data points found
        data = [data, temp];
    end
end

% plotting
b = boxchart(categorical(phase_label),data,'GroupByColor',categorical(num), 'MarkerStyle', '.');
ylabel('Facilitation')
xlabel('phase')
title('Fac of pulses per phase')
hold on

% straight lines between clusters to mark differences clearly
phase_separation = phases + (phases(2)-phases(1))/2;
xline(0.5:length(phases)+0.5);

% cleanup
clear b data phase_label num i j phase_separation

saveas(gcf,[folderpath, record.subject,num2str(record.session),'boxplot_Histogram_fac_of_each_phase.png'])
%% Compare the MEP per phase with a cluster plot
figure(11)
% find initial index of emg data
i = 1;
while savedData.emg(i).type ~= 2
    i = i+1;
end
initial_idx = i;

peaks = nan(1,length(savedData.fac_rat));
for i = 1:length(savedData.fac_rat)
    % average every two subjects
    idx = initial_idx + (i-1)*2;
    temp = mean([max(abs(savedData.emg(idx).signal(1,105:end)))  max(abs(savedData.emg(idx+1).signal(1,105:end)))]);

   
    % find peak 2 peak value
%     peaks(i) = abs(max(temp)-min(temp));
     peaks(i)=temp;
end


s = swarmchart(savedData.phases(1:length(savedData.fac_rat)), peaks, '.');
s.XJitter = 'rand';
s.XJitterWidth = 0.5;
hold on
grid on
scatter(unique(savedData.phases), MEP_s);
set(gca, 'YScale', 'log')
xticks([-2.7489   -1.9635   -1.1781   -0.3927    0.3927    1.1781    1.9635    2.7489 ])
xticklabels({'-157.5^\circ', '-112.5^\circ',  '-67.5^\circ',  '-22.5^\circ',   '22.5^\circ',   '67.5^\circ',  '112.5^\circ',  '157.5^\circ'})


legend(["Pair MEP", "Single MEP"],'Location','southoutside');
title('MEP per phase')
xlabel('phase')
ylabel('amplitude')


saveas(gcf,[folderpath, record.subject,num2str(record.session),'MEP size of_each_phase.png'])


%% PLOT THE SINGLE  and pair pusle during the traning,only avaliable for the new data  
figure(40)
[a,b]=rmmissing(savedData.temp_single_MEP) % find the singlöe pusle MEP 
c=savedData.Actions(b==0) % find the phase 

scatter(1:length(a), rmmissing(savedData.temp_single_MEP),'filled','SizeData',50)


coefficients = polyfit(1:length(a),a, 2);
yFit = polyval(coefficients, 1:length(a));

hold on;
plot(1:length(a), yFit, 'r-', 'LineWidth', 2);
grid on;


ylabel('Amplitude')
xlabel('Step')
title('Single pusle during training')
saveas(gcf,[folderpath, num2str(record.subject),num2str(record.session),'_phase_','_single pusle during the traning.png'])




figure(42)
scatter(1:length(b), peaks(1:length(b)),'filled','SizeData',50)



coefficients = polyfit(1:length(b),peaks(1:length(b)), 2);
yFit = polyval(coefficients, 1:length(b));

hold on;
plot(1:length(b), yFit, 'r-', 'LineWidth', 2);
grid on;



ylabel('Amplitude')
xlabel('Step')
title('Pair pusle during training')
saveas(gcf,[folderpath, num2str(record.subject),num2str(record.session),'_phase_','_pair pusle during the traning.png'])


figure(35)
a=rmmissing(savedData.temp_single_MEP)
s = swarmchart(c, a, 'filled','SizeData',50);
s.XJitter = 'rand';
s.XJitterWidth = 0.5;
grid on

ylabel('Amplitude')
xlabel('Phase')
title('Single pusle during training')
 set(gca, 'YScale', 'log')
    xticks(1:8)
    xticklabels({'-157.5^\circ', '-112.5^\circ',  '-67.5^\circ',  '-22.5^\circ',   '22.5^\circ',   '67.5^\circ',  '112.5^\circ',  '157.5^\circ'})
saveas(gcf,[folderpath, num2str(record.subject),num2str(record.session),'_phase_','_single pusle per phase during the traning.png'])



%% save all the parameters to excel file 
% subject information 
Expriment_date={date};
Subject_ID={savedData.subject};
Subject_Session=savedData.session;
Expriment_potocol={savedData.pusle};
Expriment_condition={savedData.condition};

%trainoptions 
Train_LR=criticOptions.LearnRate;
Train_step=trainOpts.MaxStepsPerEpisode;
Train_epsoide=trainOpts.StopTrainingValue;
Train_minibatchsize=agentOptions.MiniBatchSize;
Train_ExperienceBufferLength=agentOptions.ExperienceBufferLength;
Train_DiscountFactor=agentOptions.DiscountFactor;
Train_decay=agentOptions.EpsilonGreedyExploration.EpsilonDecay;

%stimulation paramter 
Stimulation_SMA=SMA_amplitude;
Stimulation_M1=M1_amplitude;
Fac_before_traning=fac_ratio;
MEPsss=MEP_s;
%trainprocess
% Traing_Time_perEP=str2double(extractBefore(trainStats.Information.ElapsedTime,"sec"))/(criticOptions.TrainingOpts.StopTrainingValue*60);

Newtable=table(Subject_ID,Subject_Session,Expriment_date,Expriment_potocol,Expriment_condition,Stimulation_SMA,Stimulation_M1, ...
    0,Train_LR,Train_step,Train_decay,Train_epsoide,Train_minibatchsize,Train_ExperienceBufferLength,Train_DiscountFactor,Fac_before_traning,MEPsss);

writetable(Newtable,[folderpath,savedData.subject,'_',num2str(savedData.session),'.xlsx']);





%% Random session3 after traning 
Random3

%% need a break ? 

pid=5;%change here for the best phase 
Optimize1

%% take a break 30 mins




%% record five mins resting EEG 



%% Random session4 after intervention 

Random4


%% record five mins resting EEG 


%% Optimized phase session 2
Optimize2


%% visluize the result 


save([folderpath,num2str(record.subject),'_',num2str(record.session),'workspace.mat'])



figure(33)
plotdata=[mean((savedData_R1.MEPs.paired+savedData_R2.MEPs.paired)/2) mean(savedData_R3.MEPs.paired)  mean(savedData_O1.MEPs.paired) mean(savedData_R4.MEPs.paired)  mean(savedData_O2.MEPs.paired) ]
bar(plotdata)
% xticklabels({'Random before traning','Random after traning','Optimized after traning','Random  after Post','Optimized after Post'})

xticklabels({'Before training','After training random','After training optimized','30-minutes after intervention random','30-minutes after intervention optimized'})
% ylim([0.8 1.3])
title('Paired MEP')
xlabel('sessions')
ylabel('MEPs')

saveas(gcf,[folderpath, num2str(record.subject),num2str(record.session),'_Average MEP in mutil sessions.png'])


%% only for random sessions 
save([folderpath,num2str(record.subject),'_',num2str(record.session),'workspace.mat'])

figure(55)

SPallr{1,1}=savedData_R1.MEPs.paired';
SPallr{1,2}=savedData_R2.MEPs.paired';
SPallr{1,3}=savedData_R3.MEPs.paired';
SPallr{1,4}=savedData_R4.MEPs.paired';


SPallr{2,1}=savedData_R1.MEPs.single';
SPallr{2,2}=+savedData_R2.MEPs.single';
SPallr{2,3}=savedData_R3.MEPs.single';
SPallr{2,4}=+savedData_R4.MEPs.single';


xlab={'R1','R2','R3','R4'};
col=[
102,194,165,200;
141,160,203,200];
col=col/255;

multiple_boxplot(SPallr',xlab,{'Single pulse','Pair pulse'},col')
xlabel('Phase')
ylabel('MEPs')

saveas(gcf,[folderpath, num2str(record.subject),num2str(record.session),'_Single VS paried MEPS AVG for the baseline only_all.png'])





%% plot the AVG facs 

phase_options=[-2.7489 -1.9635 -1.1781 -0.3927 0.3927 1.1781 1.9635 2.7489];


figure(18)
plotdata=[(savedData_R1.facs+savedData_R2.facs)/2 savedData_R3.facs  savedData_O1.facs savedData_R4.facs  savedData_O2.facs ]
bar(plotdata)
% xticklabels({'Random before traning','Random after traning','Optimized after traning','Random  after Post','Optimized after Post'})

xticklabels({'Before training','After training random','After training optimized','30-minutes after intervention random','30-minutes after intervention optimized'})
% ylim([0.8 1.3])
title('Faciliation ratio')
xlabel('sessions')
ylabel('Facs')


title('Average Fac of pulses per phase')

saveas(gcf,[folderpath, num2str(record.subject),num2str(record.session),'_Average fac in mutil sessions.png'])



figure(33)
plotdata=[(savedData_R1.MEPs.paired+savedData_R2.MEPs.paired)'/2; savedData_R3.MEPs.paired'; savedData_O1.MEPs.paired'; savedData_R4.MEPs.paired';  savedData_O2.MEPs.paired' ]

g1 = repmat({'Before training'},length((savedData_R1.MEPs.paired+savedData_R2.MEPs.paired)'/2),1);
g2 = repmat({'After training random'},length(savedData_R3.MEPs.paired'),1);
g3 = repmat({'After training optimized'},length(savedData_O1.MEPs.paired'),1);
g4 = repmat({'30-minutes after intervention random'},length(savedData_R4.MEPs.paired),1);
g5 = repmat({'30-minutes after intervention optimized'},length(savedData_O2.MEPs.paired),1);
g = [g1; g2; g3;g4;g5];

boxplot(plotdata,g)

title('Paired MEP')
% xlabel('sessions')
ylabel('MEPs')


saveas(gcf,[folderpath, num2str(record.subject),num2str(record.session),'_BoxplotMEP.png'])


figure(34) 
colors=rand(3,4)
meanWeight=median([(savedData_R1.MEPs.paired+savedData_R2.MEPs.paired)'/2 savedData_R3.MEPs.paired' savedData_O1.MEPs.paired' savedData_R4.MEPs.paired'  savedData_O2.MEPs.paired'])

boxchart([(savedData_R1.MEPs.paired+savedData_R2.MEPs.paired)'/2 savedData_R3.MEPs.paired' savedData_O1.MEPs.paired' savedData_R4.MEPs.paired'  savedData_O2.MEPs.paired']);
hold on
plot(meanWeight,'-o')
hold off

% xticklabels({'Random before traning','Random after traning','Optimized after traning','Random  after Post','Optimized after Post'})

xticklabels({'Before training','After training random','After training optimized','30-minutes after intervention random','30-minutes after intervention optimized'})
% ylim([0.8 1.3])
title('Paired MEP')
% xlabel('sessions')
ylabel('MEPs')

saveas(gcf,[folderpath, num2str(record.subject),num2str(record.session),'_Average MEP in mutil sessions.png'])




%% plot 
figure(35)
single_index1=savedData_R1.All_pusle(savedData_R1.All_phaseID==1)
pair_index1=savedData_R1.All_pusle(savedData_R1.All_phaseID==2)

single_index2=savedData_R2.All_pusle(savedData_R2.All_phaseID==1)
pair_index2=savedData_R2.All_pusle(savedData_R2.All_phaseID==2)

single_index3=savedData_R3.All_pusle(savedData_R3.All_phaseID==1)
pair_index3=savedData_R3.All_pusle(savedData_R3.All_phaseID==2)

single_index4=savedData_R4.All_pusle(savedData_R4.All_phaseID==1)
pair_index4=savedData_R4.All_pusle(savedData_R4.All_phaseID==2)


single_index5=savedData_O1.All_pusle(savedData_O1.All_phaseID==1)
pair_index5=savedData_O1.All_pusle(savedData_O1.All_phaseID==2)

single_index6=savedData_O1.All_pusle(savedData_O1.All_phaseID==1)
pair_index6=savedData_O1.All_pusle(savedData_O1.All_phaseID==2)


Random1all={}
Random2all={}
Random3all={}
Random4all={}

Optimized1all={}
Optimized2all={}

% fac ratio 
facsall=[]

for j=1:8

    facsall(1,j)=(savedData_R1.MEPs.paired(pair_index1==j)/savedData_R1.MEPs.single(single_index1==j)+savedData_R2.MEPs.paired(pair_index2==j)/savedData_R2.MEPs.single(single_index2==j))/2

    facsall(2,j)=savedData_R3.MEPs.paired(pair_index3==j)/savedData_R3.MEPs.single(single_index3==j)

    facsall(3,j)=savedData_R4.MEPs.paired(pair_index4==j)/savedData_R4.MEPs.single(single_index4==j)

    if j==find(phase_options==unique(savedData_O1.pairphases))

        facsall(4,j)=savedData_O1.MEPs.paired/savedData_O1.MEPs.single

        facsall(5,j)=savedData_O2.MEPs.paired/savedData_O2.MEPs.single

    else

        facsall(4,j)=0;

        facsall(5,j)=0;

    end 


end

X = categorical({'-157.5^\circ', '-112.5^\circ',  '-67.5^\circ',  '-22.5^\circ',   '22.5^\circ',   '67.5^\circ',  '112.5^\circ',  '157.5^\circ'});
X = reordercats(X,{'-157.5^\circ', '-112.5^\circ',  '-67.5^\circ',  '-22.5^\circ',   '22.5^\circ',   '67.5^\circ',  '112.5^\circ',  '157.5^\circ'});
bar(X,facsall')
% legend({'Random before training','Random after training','Optimized after training','Random  after Post','Optimized after Post'},)
legend({'Before training','After training random','Random after 30-minutes ','After training optimized','Optimized after 30-minutes'},'Location','southoutside')
% legend({'Before training','After training random','30-minutes after intervention random'},'Location','southoutside')


xlabel('Phase');
ylabel('Facilitation');

title('Fac of pulses per phase')

saveas(gcf,[folderpath, num2str(record.subject),num2str(record.session),'_fac per phase in mutil sessions.png'])
% movefile [pwd '\saved\','CHASMEC_' ,num2str(subject_number),'\','Session_' ,num2str(session_number)]   ['Z:\Experimental Data\2021-10 CHASMEC\Matlab data\',num2str(subject_number),' CHASMEC_' ,'\','Session_' ,num2str(session_number)]

%% Paired MEP all 

Mepsall={};
for j=1:8

    Mepsall{1,j}=zscore((savedData_R1.MEPs.paired(pair_index1==j)+savedData_R2.MEPs.paired(pair_index2==j))'/2)

    Mepsall{2,j}=zscore(savedData_R3.MEPs.paired(pair_index3==j)')

    Mepsall{3,j}=zscore(savedData_R4.MEPs.paired(pair_index4==j)')

    if j==find(phase_options==unique(savedData_O1.pairphases))

        Mepsall{4,j}=zscore(savedData_O1.MEPs.paired)

        Mepsall{5,j}=zscore(savedData_O2.MEPs.paired)

    end 

    if j~=find(phase_options==unique(savedData_O1.pairphases))

        Mepsall{4,j}=nan;

        Mepsall{5,j}=nan;

    end 


end


xlab={'-157.5', '-112.5',  '-67.5',  '-22.5',   '22.5',   '67.5',  '112.5',  '157.5'};
col=[228,26,28,200;
55,126,184,200;
77,175,74,200;
152,78,163,200;
255,127,0,200;];
col=col/255;

multiple_boxplot(Mepsall',xlab,{'Before training','After training random','After training optimized','Random after 30-minutes ','Optimized after 30-minutes'},col')
xlabel('Phase')
ylabel('MEPs')

saveas(gcf,[folderpath, num2str(record.subject),num2str(record.session),'_Paired MEPS each phase.png'])




%% single MEP 
Mepsall_single={};
for j=1:8

    Mepsall_single{1,j}=zscore((savedData_R1.MEPs.single(single_index1==j)+savedData_R2.MEPs.single(single_index2==j))'/2)

    Mepsall_single{2,j}=zscore(savedData_R3.MEPs.single(single_index3==j)')

    Mepsall_single{3,j}=zscore(savedData_R4.MEPs.single(single_index4==j)')

    if j==find(phase_options==unique(savedData_O1.pairphases))

        Mepsall_single{4,j}=zscore(savedData_O1.MEPs.single)

        Mepsall_single{5,j}=zscore(savedData_O2.MEPs.single)

    end 

    if j~=find(phase_options==unique(savedData_O1.pairphases))

        Mepsall_single{4,j}=nan;

        Mepsall_single{5,j}=nan;

    end 


end


xlab={'-157.5', '-112.5',  '-67.5',  '-22.5',   '22.5',   '67.5',  '112.5',  '157.5'};
col=[228,26,28,200;
55,126,184,200;
77,175,74,200;
152,78,163,200;
255,127,0,200;];
col=col/255;

multiple_boxplot(Mepsall_single',xlab,{'Before training','After training random','After training optimized','Random after 30-minutes ','Optimized after 30-minutes'},col')
xlabel('Phase')
ylabel('MEPs')

saveas(gcf,[folderpath, num2str(record.subject),num2str(record.session),'_Single MEPS each phase.png'])

% movefile 

%% AVG single and pair comparsion 
SPall{1,1}=(savedData_R1.MEPs.paired+savedData_R2.MEPs.paired)'/2
SPall{1,2}=savedData_R3.MEPs.paired'
SPall{1,3}=savedData_O1.MEPs.paired'
SPall{1,4}=savedData_R4.MEPs.paired'
SPall{1,5}=savedData_O2.MEPs.paired'


SPall{2,1}=(savedData_R1.MEPs.single+savedData_R2.MEPs.single)'/2
SPall{2,2}=savedData_R3.MEPs.single'
SPall{2,3}=savedData_O1.MEPs.single'
SPall{2,4}=savedData_R4.MEPs.single'
SPall{2,5}=savedData_O2.MEPs.single'


xlab={'R','R3','O1','R4 ','O2'};
col=[
102,194,165,200;
141,160,203,200];
col=col/255;

multiple_boxplot(SPall',xlab,{'Single pulse','Pair pulse'},col')
xlabel('Phase')
ylabel('MEPs')

saveas(gcf,[folderpath, num2str(record.subject),num2str(record.session),'_Single VS paried MEPS AVG.png'])


%% visulize the single and pair for per phase in each block 


indx=savedData_R2.All_pusle(savedData_R2.All_phaseID==1)

figure(17)
s = swarmchart(indx, savedData_R3.MEPs.paired, 'filled','LineWidth',30);
s.XJitter = 'rand';
s.XJitterWidth = 0.5;
hold on
grid on
% scatter(unique(savedData.phases), trainStats.MEP_s);
set(gca, 'YScale', 'log')
xticks(1:8)
xticklabels({'-157.5^\circ', '-112.5^\circ',  '-67.5^\circ',  '-22.5^\circ',   '22.5^\circ',   '67.5^\circ',  '112.5^\circ',  '157.5^\circ'})
title('Pair MEP random after training')
saveas(gcf,[folderpath,num2str(record.subject),'_',num2str(record.session),'EP',num2str(EPid),'_Pair MEP after trainning per phase.png'])


figure(18)
s = swarmchart(indx, savedData_R4.MEPs.paired, 'filled','LineWidth',30);
s.XJitter = 'rand';
s.XJitterWidth = 0.5;
hold on
grid on
ylim([0 10000])
% scatter(unique(savedData.phases), trainStats.MEP_s);
set(gca, 'YScale', 'log')
xticks(1:8)
xticklabels({'-157.5^\circ', '-112.5^\circ',  '-67.5^\circ',  '-22.5^\circ',   '22.5^\circ',   '67.5^\circ',  '112.5^\circ',  '157.5^\circ'})
title('Pair MEP random after 30mins')
saveas(gcf,[folderpath,num2str(record.subject),'_',num2str(record.session),'EP',num2str(EPid),'_Pair MEP random after 30 mins.png'])
close all 


figure(19)
s = swarmchart(indx, savedData_O1.MEPs.paired, 'filled','LineWidth',30);
s.XJitter = 'rand';
s.XJitterWidth = 0.5;
hold on
grid on
% scatter(unique(savedData.phases), trainStats.MEP_s);
set(gca, 'YScale', 'log')
xticks(1:8)
xticklabels({'-157.5^\circ', '-112.5^\circ',  '-67.5^\circ',  '-22.5^\circ',   '22.5^\circ',   '67.5^\circ',  '112.5^\circ',  '157.5^\circ'})
title('Pair MEP Optimized after training')
saveas(gcf,[folderpath,num2str(record.subject),'_',num2str(record.session),'EP',num2str(EPid),'_Pair MEP optimized.png'])


figure(20)
s = swarmchart(indx, savedData_O2.MEPs.paired, 'filled','LineWidth',30);
s.XJitter = 'rand';
s.XJitterWidth = 0.5;
hold on
grid on
% ylim([0 10000])
% scatter(unique(savedData.phases), trainStats.MEP_s);
set(gca, 'YScale', 'log')
xticks(1:8)
xticklabels({'-157.5^\circ', '-112.5^\circ',  '-67.5^\circ',  '-22.5^\circ',   '22.5^\circ',   '67.5^\circ',  '112.5^\circ',  '157.5^\circ'})
title('Pair MEP Optimized after 30mins')
saveas(gcf,[folderpath,num2str(record.subject),'_',num2str(record.session),'EP',num2str(EPid),'_Pair MEP optimized after 30mins.png'])


save([folderpath,num2str(record.subject),'_',num2str(record.session),'workspace.mat'])





%% save improtant data in to single mat file for later use 

%% Random 

% mkdir(['Z:\Experimental Data\2021-10 CHASMEC\Datacollection\','CHASMEC_Random\'])
 
randompath='Z:\Experimental Data\2021-10 CHASMEC\Datacollection\CHASMEC_Random\';
Datarandom.R1=savedData_R1;
Datarandom.R2=savedData_R2;
Datarandom.R3=savedData_R3;
Datarandom.R4=savedData_R4;

Datarandom.RTrain=savedData_R_training;


savename=([randompath,num2str(record.subject),'.mat']);

save(savename,'Datarandom');

%% decreased 
% mkdir(['Z:\Experimental Data\2021-10 CHASMEC\Datacollection\','CHASMEC_decreased\'])
decreasepath='Z:\Experimental Data\2021-10 CHASMEC\Datacollection\CHASMEC_decreased\';
Datadecrease.R1=savedData_R1;
Datadecrease.R2=savedData_R2;
Datadecrease.R3=savedData_R3;
Datadecrease.R4=savedData_R4;
Datadecrease.O1=savedData_O1;
Datadecrease.O2=savedData_O2;

Datadecrease.trainStats=trainStats;
Datadecrease.savedData=savedData;
Datadecrease.pid=pid;



savename=([decreasepath,num2str(record.subject),'.mat']);

save(savename,'Datadecrease');


%% enhanced data collection
% mkdir(['Z:\Experimental Data\2021-10 CHASMEC\Datacollection\','CHASMEC_enhanced\'])

%Seesion1 
enhancedpath='Z:\Experimental Data\2021-10 CHASMEC\Datacollection\CHASMEC_enhanced\';
Dataenhance.R1=savedData_R1;
Dataenhance.R2=savedData_R2;
Dataenhance.R3=savedData_R3;
Dataenhance.R4=savedData_R4;
Dataenhance.O1=savedData_O1;
Dataenhance.O2=savedData_O2;

Dataenhance.trainStats=trainStats;
Dataenhance.savedData=savedData;
Dataenhance.pid=pid;



savename=([enhancedpath,num2str(record.subject),'_',num2str(1),'.mat']);

save(savename,'Dataenhance');

% clear all

%% Session2
enhancedpath='Z:\Experimental Data\2021-10 CHASMEC\Datacollection\CHASMEC_enhanced\';
Dataenhance.R1=savedData_R1;
Dataenhance.R2=savedData_R2;
Dataenhance.R3=savedData_R3;
Dataenhance.R4=savedData_R4;
Dataenhance.O1=savedData_O1;
Dataenhance.O2=savedData_O2;

Dataenhance.trainStats=trainStats;
Dataenhance.savedData=savedData;
Dataenhance.pid=pid;



savename=([enhancedpath,num2str(record.subject),'_',num2str(2),'.mat']);

save(savename,'Dataenhance');

%% clear all
% %% plot the AVG facs notused
% figure(12)
% phase_options=[-2.7489 -1.9635 -1.1781 -0.3927 0.3927 1.1781 1.9635 2.7489];
% 
% Random1=load([folderpath,num2str(record.subject,'%03.f'),'_',num2str(record.session),'_savedData_Random1.mat'])
% Random1.savedData=Random1.savedData_R1;
% Random2=load([folderpath,num2str(record.subject,'%03.f'),'_',num2str(record.session),'_savedData_Random2.mat'])
% Random2.savedData=Random2.savedData_R2
% Random3=load([folderpath,num2str(record.subject,'%03.f'),'_',num2str(record.session),'_savedData_Random3.mat']')
% Random3.savedData=Random3.savedData_R3
% Random4=load([folderpath,num2str(record.subject,'%03.f'),'_',num2str(record.session),'_savedData_Random4.mat'])
% Random4.savedData=Random4.savedData_R4
% 
% Optimized1=load([folderpath,num2str(record.subject,'%03.f'),'_',num2str(record.session),'_savedData_optimized1.mat'])
% Optimized1.savedData=Optimized1.savedData_O1
% Optimized2=load([folderpath,num2str(record.subject,'%03.f'),'_',num2str(record.session),'_savedData_optimized2.mat'])
% Optimized2.savedData=Optimized2.savedData_O2
% 
% trainStats=load([folderpath,num2str(record.subject,'%03.f'),'_',num2str(record.session),'trainStats.mat'])
% 
% figure(18)
% plotdata=[(Random1.savedData.facs+Random2.savedData.facs)/2 Random3.savedData.facs  Optimized1.savedData.facs Random4.savedData.facs  Optimized2.savedData.facs ]
% bar(plotdata)
% % xticklabels({'Random before traning','Random after traning','Optimized after traning','Random  after Post','Optimized after Post'})
% 
% xticklabels({'Before training','After training random','After training optimized','30-minutes after intervention random','30-minutes after intervention optimized'})
% ylim([0.8 1.3])
% title('Faciliation ratio')
% xlabel('sessions')
% ylabel('Facs')
% 
% 
% title('Average Fac of pulses per phase')
% 
% saveas(gcf,[folderpath, num2str(record.subject),num2str(record.session),'_Average fac in mutil sessions.png'])
% 
% 
% 
% figure(33)
% plotdata=[mean((Random1.savedData.MEPs.paired+Random2.savedData.MEPs.paired)/2) mean(Random3.savedData.MEPs.paired)  mean(Optimized1.savedData.MEPs.paired) mean(Random4.savedData.MEPs.paired)  mean(Optimized2.savedData.MEPs.paired) ]
% bar(plotdata)
% % xticklabels({'Random before traning','Random after traning','Optimized after traning','Random  after Post','Optimized after Post'})
% 
% xticklabels({'Before training','After training random','After training optimized','30-minutes after intervention random','30-minutes after intervention optimized'})
% % ylim([0.8 1.3])
% title('Paired MEP')
% xlabel('sessions')
% ylabel('MEPs')
% 
% saveas(gcf,[folderpath, num2str(record.subject),num2str(record.session),'_Average MEP in mutil sessions.png'])
% 
% 
% 
% 
% %% plot 
% figure(13)
% single_index1=Random1.savedData.All_pusle(Random1.savedData.All_phaseID==1)
% pair_index1=Random1.savedData.All_pusle(Random1.savedData.All_phaseID==2)
% 
% single_index2=Random2.savedData.All_pusle(Random2.savedData.All_phaseID==1)
% pair_index2=Random2.savedData.All_pusle(Random2.savedData.All_phaseID==2)
% 
% single_index3=Random3.savedData.All_pusle(Random3.savedData.All_phaseID==1)
% pair_index3=Random3.savedData.All_pusle(Random3.savedData.All_phaseID==2)
% 
% single_index4=Random4.savedData.All_pusle(Random4.savedData.All_phaseID==1)
% pair_index4=Random4.savedData.All_pusle(Random4.savedData.All_phaseID==2)
% 
% Random1all={}
% Random2all={}
% Random3all={}
% Random4all={}
% 
% Optimized1all={}
% Optimized2all={}
% 
% facsall=[]
% 
%      for j=1:8
% 
%      facsall(1,j)=(Random1.savedData.MEPs.paired(pair_index1==j)/Random1.savedData.MEPs.single(single_index1==j)+Random2.savedData.MEPs.paired(pair_index2==j)/Random2.savedData.MEPs.single(single_index2==j))/2  
%    
% %      facsall(2,j)=
%   
%      facsall(2,j)=Random3.savedData.MEPs.paired(pair_index3==j)/Random3.savedData.MEPs.single(single_index3==j)  
%      
%      facsall(3,j)=Random4.savedData.MEPs.paired(pair_index4==j)/Random4.savedData.MEPs.single(single_index4==j)  
%   
% %      facsall(4,j)=Optimized1.savedData.MEPs.paired(pair_index==j)/Optimized1.savedData.MEPs.single(single_index==j)  
% %    
% %      facsall(5,j)=Optimized2.savedData.MEPs.paired(pair_index==j)/Optimized2.savedData.MEPs.single(single_index==j)  
%    
%      end 
% 
% 
% 
% 
% X = categorical({'-157.5^\circ', '-112.5^\circ',  '-67.5^\circ',  '-22.5^\circ',   '22.5^\circ',   '67.5^\circ',  '112.5^\circ',  '157.5^\circ'});
% X = reordercats(X,{'-157.5^\circ', '-112.5^\circ',  '-67.5^\circ',  '-22.5^\circ',   '22.5^\circ',   '67.5^\circ',  '112.5^\circ',  '157.5^\circ'});
% bar(X,facsall')
% % legend({'Random before training','Random after training','Optimized after training','Random  after Post','Optimized after Post'},)
% % legend({'Before training','After training random','After training optimized','30-minutes after intervention random','30-minutes after intervention optimized'},'Location','southoutside')
% legend({'Before training','After training random','30-minutes after intervention random'},'Location','southoutside')
% 
% 
% xlabel('Phase');
% ylabel('Facilitation');
% 
% title('Fac of pulses per phase')
% 
% saveas(gcf,[folderpath, num2str(record.subject),num2str(record.session),'_fac per phase in mutil sessions.png'])
% 
% % close all
% 
% titletex={'Before training','After training random','After training optimized','30-minutes after intervention random','30-minutes after intervention optimized'}
% 
% 
% 
% 
% %% organize the data for mean phase
% facs_all={}
% facs_final_mean={};
% facs_all_std={};
% data_facs_all={Random1,Random2,Random3,Optimized1,Random4,Optimized2}
% for i=1:6
%   
%      data_temp=data_facs_all{i}
%      for j=1:8
% 
%      facs_all{j,1}=mean(data_temp.savedData.MEPs.paired(pair_index==j))
%      facs_all{j,2}=mean(data_temp.savedData.MEPs.single(single_index==j)) 
% 
%      facs_all_s{j,1}=std(data_temp.savedData.MEPs.paired(pair_index==j))
%      facs_all_s{j,2}=std(data_temp.savedData.MEPs.single(single_index==j)) 
%    
%      end 
%     facs_final_mean{i}=facs_all;
%     facs_all_std{i}=facs_all_s;
% 
% end    
% 
% % Plot the average of each phase for single and pair MEP 
% facs_final_means={};
% Randomcom_single=(cell2mat(facs_final_mean{1}(:,1))+cell2mat(facs_final_mean{2}(:,1)))/2
% Randomcom_pair=(cell2mat(facs_final_mean{1}(:,2))+cell2mat(facs_final_mean{2}(:,2)))/2
% facs_final_means{1}(:,1)=num2cell(Randomcom_single)
% facs_final_means{1}(:,2)=num2cell(Randomcom_pair)
% facs_final_means{2}=facs_final_mean{3};
% facs_final_means{3}=facs_final_mean{4};
% facs_final_means{4}=facs_final_mean{5};
% facs_final_means{5}=facs_final_mean{6};
% 
% 
% 
% 
% 
% for m=1:5
% single_index1=[]
% pair_index1=[]
% single_mep=[]
% pair_mep=[]
%     dataphase=facs_final_mean{m};
%   
% 
% 
% 
%     for t=1:8
%         single_index1=[single_index1,t*ones(1,length(dataphase{t,1}))]
% 
%         single_mep=[single_mep,dataphase{t,1}]
% 
%         pair_index1=[pair_index1,t*ones(1,length(dataphase{t,2}))]
% 
%         pair_mep=[pair_mep,dataphase{t,2}]
% 
%     end
% 
%     close all;
%     s = swarmchart(single_index1, single_mep, 'filled','SizeData',100);
%     s.XJitter = 'rand';
%     s.XJitterWidth = 0.5;
%     hold on
%     s = swarmchart(pair_index1, pair_mep, 'filled','SizeData',100);
%     hold on
%     grid on
% 
%     scatter(1:8, trainStats.trainStats.MEP_s,'filled','SizeData',100, 'MarkerEdgeColor','flat');
% 
%     set(gca, 'YScale', 'log')
%     xticks(1:8)
%     xticklabels({'-157.5^\circ', '-112.5^\circ',  '-67.5^\circ',  '-22.5^\circ',   '22.5^\circ',   '67.5^\circ',  '112.5^\circ',  '157.5^\circ'})
% 
% 
%     legend(["Single MEP", "Pair MEP",'Single\_Avg'],'Location','southoutside');
%     title(titletex{m})
%     xlabel('phase')
%     ylabel('amplitude')
%  set(gcf,'position',[200,300,800,600]);
%     saveas(gcf,[folderpath, num2str(record.subject),num2str(record.session),'_Session_',num2str(m),'_Mean_fac per phase in one seesions.png'])
% 
%       
%            close all;
% 
% end
% 
% 
% 
% 
% 
% 
% %% organize the data for each phase 
% facs_all={}
% facs_final={}
% 
% % data_facs_all={Random1,Random2,Random3,Optimized1,Random4,Optimized2}
% data_facs_all={Random1,Random2,Random3,Optimized1,Random4,Optimized2}
% for i=1:6
%   
%      data_temp=data_facs_all{i}
%      for j=1:8
% 
%      facs_all{j,1}=data_temp.savedData.MEPs.paired(pair_index==j)
%      facs_all{j,2}=data_temp.savedData.MEPs.single(single_index==j)  
%    
%      end 
%     facs_final{i}=facs_all;
% end    
% 
% 
% % plot the mep per phase for each session 
% titletex={'Before training1','Before training2','After training random','After training optimized','30-minutes after intervention random','30-minutes after intervention optimized'}
% 
% 
% for m=1:6
%     single_index1=[]
% pair_index1=[]
% single_mep=[]
% pair_mep=[]
%     dataphase=facs_final{m};
% 
%     for t=1:8
%         single_index1=[single_index1,t*ones(1,length(dataphase{t,1}))]
% 
%         single_mep=[single_mep,dataphase{t,1}]
% 
%         pair_index1=[pair_index1,t*ones(1,length(dataphase{t,2}))]
% 
%         pair_mep=[pair_mep,dataphase{t,2}]
% 
%     end
% 
%     close all;
%     s = swarmchart(single_index1, single_mep, 'filled','SizeData',100);
%     s.XJitter = 'rand';
%     s.XJitterWidth = 0.5;
%     hold on
%     s = swarmchart(pair_index1, pair_mep, 'filled','SizeData',100);
%     hold on
%     grid on
% 
%     scatter(1:8, trainStats.trainStats.MEP_s,'filled','SizeData',100, 'MarkerEdgeColor','flat');
% 
%     set(gca, 'YScale', 'log')
%     xticks(1:8)
%     xticklabels({'-157.5^\circ', '-112.5^\circ',  '-67.5^\circ',  '-22.5^\circ',   '22.5^\circ',   '67.5^\circ',  '112.5^\circ',  '157.5^\circ'})
% 
% 
%     legend(["Single MEP", "Pair MEP",'Single\_Avg'],'Location','southoutside');
%     title(titletex{m})
%     xlabel('phase')
%     ylabel('amplitude')
%  set(gcf,'position',[200,300,800,600]);
%     saveas(gcf,[folderpath, num2str(record.subject),num2str(record.session),'_Session_',num2str(m),'_fac per phase in one seesions.png'])
% 
%       
%          
% 
% end
% 
% 
% %% plot the mean mep per phase for all session tegether 
% titletex={'Before training','After training random','After training optimized','30-minutes after intervention random','30-minutes after intervention optimized'}
% 
% for t=1:8
% 
% single_index1=[]
% pair_index1=[]
% single_mep=[]
% pair_mep=[]
% 
% for m=1:5
%     dataphase=facs_final_means{m};
%         single_index1=[single_index1,m*ones(1,length(dataphase{t,1}))]
% 
%         single_mep=[single_mep,dataphase{t,1}]
% 
%         pair_index1=[pair_index1,m*ones(1,length(dataphase{t,2}))]
% 
%         pair_mep=[pair_mep,dataphase{t,2}]
% 
% end
% 
% 
%      close all;
%     s = swarmchart(single_index1, single_mep, 'filled','SizeData',100);
%     s.XJitter = 'rand';
%     s.XJitterWidth = 0.5;
%     hold on
%     s = swarmchart(pair_index1, pair_mep, 'filled','SizeData',100);
%     hold on
%     grid on
% 
% %     scatter(t, trainStats.trainStats.MEP_s(t),'filled','SizeData',100, 'MarkerEdgeColor','flat');
% %     yline(trainStats.trainStats.MEP_s(t))
% 
%     yline(trainStats.trainStats.MEP_s(t),'-','Initial single MEP');
% 
%     set(gca, 'YScale', 'log')
%     xticks(1:6)
%     xticklabels(titletex)
% 
% 
%     legend(["Single MEP", "Pair MEP",'Initial single MEP'],'Location','southoutside');
%     title(num2str(rad2deg(phase_options(t))))
%     xlabel('phase')
%     ylabel('Amplitude')
% 
%    set(gcf,'position',[200,300,800,600]);
% 
%     saveas(gcf,[folderpath num2str(record.subject),num2str(record.session),'_phase_',num2str(t),'_fac only in one in mutil sessions.png'])
%     close all;
% 
% 
% 
% end 
% 
% 
% %% plot the mep per phase for all session tegether 
% 
% 
% for t=1:8
% 
% single_index1=[]
% pair_index1=[]
% single_mep=[]
% pair_mep=[]
% 
% for m=1:6
%     dataphase=facs_final{m};
%         single_index1=[single_index1,m*ones(1,length(dataphase{t,1}))]
% 
%         single_mep=[single_mep,dataphase{t,1}]
% 
%         pair_index1=[pair_index1,m*ones(1,length(dataphase{t,2}))]
% 
%         pair_mep=[pair_mep,dataphase{t,2}]
% 
% end
% 
% 
%      close all;
%     s = swarmchart(single_index1, single_mep, 'filled','SizeData',100);
%     s.XJitter = 'rand';
%     s.XJitterWidth = 0.5;
%     hold on
%     s = swarmchart(pair_index1, pair_mep, 'filled','SizeData',100);
%     hold on
%     grid on
% 
%     scatter(t, trainStats.trainStats.MEP_s(t),'filled','SizeData',100, 'MarkerEdgeColor','flat');
% 
%     set(gca, 'YScale', 'log')
%     xticks(1:6)
%     xticklabels(titletex)
% 
% 
%     legend(["Single MEP", "Pair MEP",'Single\_Avg'],'Location','southoutside');
%     title(num2str(rad2deg(phase_options(t))))
%     xlabel('phase')
%     ylabel('amplitude')
% 
%    set(gcf,'position',[200,300,800,600]);
% 
%     saveas(gcf,[folderpath, num2str(record.subject),num2str(record.session),'_phase_',num2str(t),'_all points_fac only in one in mutil sessions.png'])
%     close all;
% 
% 
% 
% end 
% 
% 
% 
% 
% 
% 
% 
% 
% 
