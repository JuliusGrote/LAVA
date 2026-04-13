import random
import numpy as np

random.seed(42)

n_trials = 600

phase_conditions = [0, np.pi, np.nan]

data_save = "C:/Users/Julius/Meine Ablage/Studium/GTC/LabRotations/ZiemannLab/LAVA/data/"

# A multiplier of 2 means each block has 6 trials (2 zero, 2 trough, 2 random).
multiplier = 2 
block = phase_conditions * multiplier 
num_blocks = int(n_trials / len(block))

experiment_sequence = []

# Generate the sequence block by block
for _ in range(num_blocks):
    current_block = block.copy()
    random.shuffle(current_block)
    experiment_sequence.extend(current_block)

# exchange np.nan with random phase values between 0 and 2*pi

experiment_sequence = [random.uniform(0, 2 * np.pi) if np.isnan(phase) else phase for phase in experiment_sequence]

with open(f"{data_save}phase_schedule.csv", "w") as f:
    for values in experiment_sequence:
        f.write(f"{values},\n")