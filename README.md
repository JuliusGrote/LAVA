# LAVA (LAbjack VAlidation) Project
The project aimes to compare the extrablished bossdevice closed-loop setup with the Neurosimo loop.

## Directories
The 'code' folder hosts all scripts relevant for running the actual experiment. 'Data analysis' conversly hosts the matlab and python code to calculate and plot the statistical results.
It uses 'data' where all collected data as well as the random phases and filter coefficients are stored which are used by the scripts in 'code'.
The graphics produced can be found in the 'figures' folder. All relevant plugins/toolboxes that deliver external functions are in 'toolboxes'. Statistical analysis results are stored in 'stats_results'.

### Control PC location
The LAVA on the control PC can be found on the desktop/LAVA folder. It only contains the bossdevice toolbox for the MATLAB API and UDP listener script. 

### Neurosimo PC location
Its counterpart is found in the home/projects/LAVA directory. The folder is stuctured as a NeuroSimo compatible project. 

## Running the dual closed loop
1. Ensure that this folder exists on both the Control PC and Neurosimo PC. 

2. Turn on NeurOne Black, battery and amplifier. Start the EEG measurement and recording in the NeurOne software on the Control  PC (use the 20205-10LOOP protocol).

3. Turn on the BOSS device. Start Matlab and run the LAVA_BOSS.m script. After a certain while/section the UDP listener will be started. This means that the script is now waiting for the start signal from the Neurosimo PC.

4. On the Neurosimo PC start Neurosimo (use the existing LAVA project with the LAVA_Neurosimo.py   decider and the preprocessor turned off). Click start session. Now the experiment starts. 
In the NeurOne application you should now see triggers from the BOSS ('10', '20', '30') and Neurosimo ('Stimulation Out'). 

5. Wait until all triggers are given (by both loops!) and then end the recording and measurement and, if appropriate, shut down all scripts and devices. 

## Adjustments and function
If you want to adjust certain parameters do so in the beginning of the python Neurosimo decider and in the respective section of the matlab script (properties of the 'bd' object). 
Both the loops are started at the same time by a udp trigger sent from the Neurosimo PC to the Control PC.
The source code can be found in the process() of the LAVA_Neurosimo.py file and the UDP listener is implmemented in 'Matlab_udp_listener'.
