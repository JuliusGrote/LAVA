"""
Phastimate decider module for NeuroSimo.

Reference:
Zrenner, C., Galevska, D., Nieminen, J.O., Baur, D., Stefanou, M.I., Ziemann, U. (2020). 
The shaky ground truth of real-time phase estimation. NeuroImage, 214, 116761.
https://doi.org/10.1016/j.neuroimage.2020.116761

Available at:
https://www.sciencedirect.com/science/article/pii/S1053811920302262

MATLAB implementation:
https://github.com/bnplab/phastimate

This version is based on MATLAB adaptation of Phastimate by Joonas Laurinoja at Aalto University.
"""

import csv
import subprocess
from typing import Dict, List, Optional, Tuple, Union
import numpy as np
from scipy.signal import filtfilt, hilbert
from scipy.io import loadmat
from spectrum import aryule

SUBJECT_ID = 'Sub-01' # increase iteravely for each subject


# Load MATLAB coefficients directly as numpy array
mat_data = loadmat('filter_coeffs.mat') # correct, matches BOSSDevice
BANDPASS_FILTER_COEFFICIENTS = np.array(mat_data['coeffs'].flatten()) # correct, matches BOSSDevice


# EEG channel indices for C3 referencing matches BOSSDevice
C3_CHANNEL_INDEX = 4 # add 1 to convert to MATLAB (1-based)
REFERENCE_CHANNEL_INDICES = [20, 22, 24, 26]  # Reference channels for C3, again 0-based
REFERENCE_WEIGHT = 0.25 

# Phase estimation constants
TARGET_PHASE_RADIANS = [np.pi, 0, 'random']  # Three phase targets: Trough (π), Peak (0), Random (from CSV)
RANDOM_PHASES_CSV_PATH = 'data/random_phases.csv'  # Path to CSV file with random phase targets
DEFAULT_PHASE_TOLERANCE = np.pi / 40  # Single tolerance value for all phases

# Phastimate algorithm parameters (Zrenner et al. 2020)
# Note: BOSSdevice firmware parameters are compiled in Simulink model (not accessible via API)
# BOSSdevice examples show offset_samples=3-6 (depends on loop delay)
# Zrenner2020: Uses 500ms windows, 128-sample Hilbert window at 500Hz (after 10x decimation from 5kHz)
DEFAULT_HILBERT_WINDOW_SIZE = 128  
DEFAULT_EDGE_SAMPLES = 64
# Note: BOSS offset_samples=3-6 at 500Hz, scaled to your sampling rate this equals ~30-60 at 5000Hz pre-decimation
DEFAULT_AR_MODEL_ORDER = 25  # AR model order - good balance between under/overfitting
DEFAULT_DOWNSAMPLE_RATIO = 10  # Matches BOSSDevice: 5000Hz -> 500Hz

# Processing timing constants
DEFAULT_PROCESSING_INTERVAL_SECONDS = 0.05

DEFAULT_BUFFER_SIZE_SECONDS = 0.5  # ?? 
TRIGGER_COOLDOWN_SECONDS = 2.0  # Minimum time between triggers in s - correct (set here, can be adjusted for BOSS)

# Set the Number of Trials per phase
N_TRIALS_PER_PHASE = 200  # 200 trials per phase
N_PHASES = 3  # Number of phase targets (π, 0, random)

class Decider:
    """
    Real-time EEG phase estimation and trigger scheduling using Phastimate algorithm.
    
    This class implements the Phastimate algorithm for real-time phase estimation
    of EEG signals and schedules triggers based on target phase detection.
    """
    
    def __init__(self, num_eeg_channels: int, num_emg_channels: int, sampling_frequency: float):
        """
        Initialize the Decider with parameters and filter design.
        
        Args:
            num_eeg_channels: Number of EEG channels (unused but kept for interface compatibility)
            num_emg_channels: Number of EMG channels (unused but kept for interface compatibility)
            sampling_frequency: Sampling frequency in Hz
        """

        self.subject_id = SUBJECT_ID
        self.sampling_frequency = sampling_frequency

        # Phastimate algorithm parameters
        self.hilbert_window_size = DEFAULT_HILBERT_WINDOW_SIZE
        self.edge_samples = DEFAULT_EDGE_SAMPLES
        self.ar_model_order = DEFAULT_AR_MODEL_ORDER
        self.downsample_ratio = DEFAULT_DOWNSAMPLE_RATIO

        # Processing timing parameters
        self.processing_interval_seconds = DEFAULT_PROCESSING_INTERVAL_SECONDS
        self.buffer_size_seconds = DEFAULT_BUFFER_SIZE_SECONDS

        # Convert timing parameters to samples
        self.processing_interval_samples = int(self.processing_interval_seconds * sampling_frequency)
        self.buffer_size_samples = int(self.buffer_size_seconds * sampling_frequency)
        print(f"Buffer size in samples: {self.buffer_size_samples}")


        # Filter coefficients
        self.bandpass_filter_coefficients = BANDPASS_FILTER_COEFFICIENTS

        # Phase targeting parameters
        self.target_phases = TARGET_PHASE_RADIANS  # List of all phase targets
        self.phase_tolerance = DEFAULT_PHASE_TOLERANCE
        
        
        
        # Maximum number of future samples to consider for trigger scheduling
        self.max_future_samples = int(self.edge_samples / 2)

        # Leave empty when mTMS device is not used
        self.targets = []
        

        # Initialize logs - single structure for all phases
        self.timestamp_log = []
        self.triggertimes_log = []
        self.phases_log = []  # Track which phase was targeted

        # Number of Trials per phase
        self.n_trials_per_phase = N_TRIALS_PER_PHASE
        self.n_phases = N_PHASES
        # Number of warm-up rounds to prevent first-call delays (see README.md for details)
        self.warm_up_rounds = 2

        # Phase tracking
        self.current_phase_index = 0  # Track which phase we're currently on
        self.current_phase_trials = 0  # Track trials completed in current phase
        
        # Trigger cooldown tracking (2 seconds minimum between triggers)
        self.last_trigger_time = None
        self.trigger_cooldown_seconds = TRIGGER_COOLDOWN_SECONDS
        
        # Load random phases from CSV file
        self.random_phases = self._load_random_phases(RANDOM_PHASES_CSV_PATH)
        self.random_phase_index = 0  # Track current position in random phases list


        # Set initial target phase (handle random case)
        if self.target_phases[0] == 'random':
            self.target_phase_radians = self.random_phases[self.random_phase_index]
        else:
            self.target_phase_radians = self.target_phases[0]
        
        # Total trials across all phases
        total_trials = self.n_trials_per_phase * self.n_phases
        self.trials_left = total_trials


    def get_configuration(self) -> Dict[str, Union[int, bool, List]]:
        """
        Return the configuration for the processing interval and sample window.
        
        Returns:
            Dictionary containing processing configuration parameters
        """
        # Call MATLAB script on remote Windows PC
        try:
            subprocess.run([
                'ssh', 'BNP Brown@192.168.1.100',  # Replace with actual IP
                'matlab.exe', '-batch', '"run(\'C:\\path\\to\\script.m\'); exit;"'
            ], capture_output=True, timeout=30)
        except Exception as e:
            print(f"Warning: Remote MATLAB call failed: {e}")
        
        return {
            'processing_interval_in_samples': self.processing_interval_samples,
            'process_on_trigger': False,
            'sample_window': [-(self.buffer_size_samples - 1), 0],
            'events': [],
            'sensory_stimuli': [],
        }

    def _load_random_phases(self, csv_path: str) -> np.ndarray:
        """
        Load random phase targets from a CSV file.
        
        Args:
            csv_path: Path to the CSV file containing random phase values (in radians)
            
        Returns:
            Array of random phase values
        """
        phases_array = None
        generate_phases = False
        try:
            with open(csv_path, 'r') as f:
                reader = csv.reader(f)

                phases = [float(row[0]) for row in reader]
                
            if not phases:
                print(f"Warning: No phases found in {csv_path}. Generating random phases.")
                # Generate default random phases if file is empty
                generate_phases = True
            else:
                phases_array = np.array(phases)
                print(f"Loaded {len(phases_array)} random phase targets from {csv_path}")
                
                # Ensure we have enough phases for all trials
            if len(phases_array) < self.n_trials_per_phase:
                print(f"Warning: CSV has {len(phases_array)} phases but need {self.n_trials_per_phase}. Generating random phases.")
                # Scrub old phases and generate new ones
                generate_phases = True

        except FileNotFoundError:
            print(f"Warning: Random phases CSV file not found at {csv_path}. Generating random phases.")
            generate_phases = True
        except Exception as e:
            print(f"Error loading random phases from CSV: {e}. Generating random phases.")
            generate_phases = True
        
        # Always save the phases array to ensure newly generated phases are saved
        if generate_phases:
            np.random.seed(42)  # For reproducibility
            phases_array = np.random.uniform(-np.pi, np.pi, self.n_trials_per_phase)
            print(f"Generated {len(phases_array)} random phase targets.")
            save_path = csv_path
            with open(save_path, 'w', newline='') as f:
                writer = csv.writer(f)
                for phase in phases_array:
                    writer.writerow([phase])
            print(f"Saved {len(phases_array)} phases to {save_path}")
        
        return phases_array

    def process(self, current_time: float, timestamps: np.ndarray, valid_samples: np.ndarray, 
                eeg_buffer: np.ndarray, emg_buffer: np.ndarray, current_sample_index: int, 
                ready_for_trial: bool, is_trigger: bool, is_event: bool, event_type: str, is_coil_at_target: bool) -> Optional[Dict[str, float]]:
        """
        Process the EEG data to estimate phase and schedule a trigger.
        
        Args:
            current_time: Current timestamp
            timestamps: Array of sample timestamps
            valid_samples: Boolean array indicating valid samples
            eeg_buffer: EEG data buffer (samples x channels)
            emg_buffer: EMG data buffer (unused)
            current_sample_index: Current sample index
            ready_for_trial: Whether the system is ready for a new trial
            is_trigger: Whether a trigger event occurred
            is_event: Whether an event occurred
            event_type: Type of event
            
        Returns:
            Dictionary with 'timed_trigger' key and execution time, or None if no trigger scheduled
        """
        # Check if we need to advance to next phase
        if self.current_phase_trials >= self.n_trials_per_phase:
            self.current_phase_index += 1
            self.current_phase_trials = 0
            
            if self.current_phase_index >= self.n_phases:
                # All phases completed
                print("\n=== All phases completed! ===")
                self._save_logs(self.subject_id)

                # Stop the session using ROS 2 service
                print("\n---------------------------- Attempting to stop session via ROS service... ----------------------------\n")
                try:
                    result = subprocess.run(
                        ["ros2", "service", "call", "/system/session/stop", "system_interfaces/srv/StopSession"],
                        capture_output=True,
                        text=True,
                        timeout=5
                    )
                    if result.returncode == 0:
                        print("Session stopped successfully!")
                    else:
                        print(f"Failed to stop session: {result.stderr}")
                except Exception as e:
                    print(f"Error stopping session: {e}")
                return None
            
            # Update target phase for new phase
            phase_names = ["Trough (π)", "Peak (0)", "Random (from CSV)"]
            print(f"\n=== Starting Phase {self.current_phase_index + 1}: {phase_names[self.current_phase_index]} ===")
            
            # Set target phase (handle random case)
            if self.target_phases[self.current_phase_index] == 'random':
                # For random phase, reset the index when starting random phase block
                self.random_phase_index = 0
                self.target_phase_radians = self.random_phases[self.random_phase_index]
            else:
                self.target_phase_radians = self.target_phases[self.current_phase_index]
        
        # For random phase, use the next value from the loaded CSV
        if self.target_phases[self.current_phase_index] == 'random':
            self.target_phase_radians = self.random_phases[self.random_phase_index]
            self.random_phase_index += 1
            print(f"Random phase target from CSV: {self.target_phase_radians:.4f} rad ({np.degrees(self.target_phase_radians):.2f}°)")
        
        # Check trigger cooldown - skip trial attempts during cooldown
        if self.last_trigger_time is not None:
            time_since_last_trigger = current_time - self.last_trigger_time
            if time_since_last_trigger < self.trigger_cooldown_seconds:
                # Don't even attempt processing during cooldown period
                return None
        
        # Early returns for invalid states
        if not ready_for_trial or not np.all(valid_samples):
            return None

        # Extract C3 channel with common average reference
        c3_referenced_data = self._extract_c3_referenced_data(eeg_buffer)
        if c3_referenced_data is None:
            return None

        # Preprocess the data
        preprocessed_data = self._preprocess_eeg_data(c3_referenced_data)

        # Estimate future phases using Phastimate algorithm
        estimated_phases = self._estimate_phases(preprocessed_data)
        if estimated_phases is None:
            return None

        # Find optimal trigger timing
        trigger_timing = self._find_optimal_trigger_timing(estimated_phases, current_time)
        
        # Check trigger cooldown - prevent triggers within 2 seconds of last trigger
        if trigger_timing is not None and self.last_trigger_time is not None:
            time_between_triggers = trigger_timing['timed_trigger'] - self.last_trigger_time
            if time_between_triggers < self.trigger_cooldown_seconds:
                print(f"Cooldown violation: {time_between_triggers:.3f}s < {self.trigger_cooldown_seconds:.1f}s - trigger rejected")
                trigger_timing = None
        
        # Update last trigger time if a trigger was scheduled
        if trigger_timing is not None:
            self.last_trigger_time = trigger_timing['timed_trigger']  # Use actual trigger execution time

        # log timestamp, trigger time, and phase
        self.timestamp_log.append(current_time)
        self.triggertimes_log.append(trigger_timing)
        self.phases_log.append(self.target_phase_radians)
        
        # Update counters
        self.current_phase_trials += 1
        self.trials_left -= 1
        
        phase_names = ["Trough (π)", "Peak (0)", "Random (from CSV)"]
        print(f"Phase {self.current_phase_index + 1} - Trial {self.current_phase_trials} "
              f"(global: {self.n_trials_per_phase * self.n_phases - self.trials_left + 1})")

        return trigger_timing

    def _extract_c3_referenced_data(self, eeg_buffer: np.ndarray) -> Optional[np.ndarray]:
        """
        Extract C3 channel data with common average reference.
        
        Args:
            eeg_buffer: EEG data buffer (samples x channels)
            
        Returns:
            Referenced C3 data or None if extraction fails
        """
        try:
            c3_data = eeg_buffer[:, C3_CHANNEL_INDEX]
            reference_data = np.sum(eeg_buffer[:, REFERENCE_CHANNEL_INDICES], axis=1)
            return c3_data - REFERENCE_WEIGHT * reference_data
        except IndexError:
            print("Error: EEG buffer does not have expected number of channels")
            return None

    def _preprocess_eeg_data(self, data: np.ndarray) -> np.ndarray:
        """
        Preprocess EEG data by demeaning and downsampling.
        
        Args:
            data: Input EEG data
            
        Returns:
            Preprocessed and downsampled data
        """
        # Remove DC component
        demeaned_data = data - np.mean(data)
        
        # Downsample the data
        return demeaned_data[::self.downsample_ratio]

    def _estimate_phases(self, data: np.ndarray) -> Optional[np.ndarray]:
        """
        Estimate future phases using the Phastimate algorithm.
        
        Args:
            data: Preprocessed EEG data
            
        Returns:
            Array of estimated phases or None if estimation fails
        """
        estimated_phases, _ = self.phastimate(
            data,
            self.bandpass_filter_coefficients, 
            [1.0], 
            self.edge_samples, 
            self.ar_model_order, 
            self.hilbert_window_size
        )
        
        return estimated_phases

    def _find_optimal_trigger_timing(self, estimated_phases: np.ndarray, current_time: float) -> Optional[Dict[str, float]]:
        """
        Find optimal trigger timing based on estimated phases.
        
        Args:
            estimated_phases: Array of estimated phase values
            current_time: Current timestamp
            
        Returns:
            Dictionary with trigger timing or None if no suitable timing found
        """
        # Extract future phase estimates (second half of the estimation window)
        num_samples = estimated_phases.shape[0]
        future_phase_estimates = estimated_phases[num_samples // 2:]
        
        # Limit to maximum future samples
        future_phase_estimates = future_phase_estimates[:self.max_future_samples]

        # Calculate phase differences from target
        phase_differences = np.angle(np.exp(1j * (future_phase_estimates - self.target_phase_radians)))

        # Find the sample with minimum phase difference
        optimal_sample_index = np.argmin(np.abs(phase_differences))
        min_phase_difference = phase_differences[optimal_sample_index]

        # Check if phase difference is within tolerance
        if np.abs(min_phase_difference) > self.phase_tolerance:
            print(f'Phase difference exceeds tolerance: {min_phase_difference:.3f} radians')
            return None

        # Calculate trigger execution time
        time_offset_seconds = (optimal_sample_index * self.downsample_ratio) / self.sampling_frequency
        execution_time = current_time + time_offset_seconds

        print(f'Trigger scheduled {time_offset_seconds:.3f} seconds from now')

        return {'timed_trigger': execution_time}

    def phastimate(self, data: np.ndarray, filter_b: np.ndarray, filter_a: List[float], 
                   edge_samples: int, ar_order: int, hilbert_window_size: int,
                   offset_correction: int = 0, iterations: Optional[int] = None, 
                   ar_method: str = 'aryule') -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
        """
        Estimate the phase of the EEG signal using autoregressive modeling and Hilbert transform.
        
        This is the core Phastimate algorithm that performs:
        1. Bandpass filtering of the input signal
        2. Autoregressive (AR) modeling for forward prediction
        3. Hilbert transform for phase extraction
        
        Args:
            data: Input EEG signal
            filter_b: Numerator coefficients of the bandpass filter
            filter_a: Denominator coefficients of the bandpass filter  
            edge_samples: Number of edge samples to remove after filtering
            ar_order: Order of the autoregressive model
            hilbert_window_size: Size of the window for Hilbert transform
            offset_correction: Offset correction (unused)
            iterations: Number of forward prediction iterations
            ar_method: Method for AR parameter estimation ('aryule')
            
        Returns:
            Tuple of (estimated_phases, estimated_amplitudes) or (None, None) if estimation fails
        """
        # Calculate default number of iterations if not specified
        if iterations is None:
            iterations = edge_samples + int(np.ceil(hilbert_window_size / 2))

        # Validate data length for filtering
        min_padding_length = 3 * (max(len(filter_a), len(filter_b)) - 1)
        if data.shape[0] <= min_padding_length:
            print(f"Insufficient data for filtering: {data.shape[0]} <= {min_padding_length}")
            return None, None

        # Apply bandpass filter
        filtered_data = filtfilt(filter_b, filter_a, data)

        # Remove edge samples to mitigate filter transients
        if filtered_data.shape[0] <= 2 * edge_samples:
            print(f"Insufficient data after edge removal: {filtered_data.shape[0]} <= {2 * edge_samples}")
            return None, None

        edge_removed_data = filtered_data[edge_samples:-edge_samples]

        # Fit autoregressive model
        ar_coefficients = self._fit_ar_model(edge_removed_data, ar_order, ar_method)
        if ar_coefficients is None:
            return None, None

        # Perform forward prediction
        predicted_data = self._forward_predict(edge_removed_data, ar_coefficients, iterations)

        # Extract phase and amplitude using Hilbert transform
        return self._extract_phase_amplitude(predicted_data, hilbert_window_size)

    def _fit_ar_model(self, data: np.ndarray, ar_order: int, ar_method: str) -> Optional[np.ndarray]:
        """
        Fit autoregressive model to the data.
        
        Args:
            data: Input data for AR modeling
            ar_order: Order of the AR model
            ar_method: Method for AR parameter estimation
            
        Returns:
            AR coefficients or None if fitting fails
        """
        if len(data) < ar_order:
            print(f"Insufficient data for AR model: {len(data)} < {ar_order}")
            return None

        if ar_method == 'aryule':
            try:
                ar_params, _, _ = aryule(data, ar_order)
                # Flip and negate coefficients for prediction equation
                return -1 * ar_params[::-1]
            except Exception as e:
                print(f"AR model fitting failed: {e}")
                return None
        else:
            raise ValueError(f'Unknown AR method: {ar_method}')

    def _forward_predict(self, data: np.ndarray, ar_coefficients: np.ndarray, iterations: int) -> np.ndarray:
        """
        Perform forward prediction using AR model.
        
        Args:
            data: Historical data
            ar_coefficients: AR model coefficients
            iterations: Number of prediction steps
            
        Returns:
            Extended data with predictions
        """
        # Initialize prediction array
        total_length = len(data) + iterations
        predicted_data = np.zeros(total_length)
        predicted_data[:len(data)] = data

        # Perform iterative prediction
        ar_order = len(ar_coefficients)
        for i in range(iterations):
            prediction_index = len(data) + i
            data_window = predicted_data[prediction_index - ar_order:prediction_index]
            predicted_data[prediction_index] = np.sum(ar_coefficients * data_window)

        return predicted_data

    def _extract_phase_amplitude(self, data: np.ndarray, window_size: int) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
        """
        Extract phase and amplitude using Hilbert transform.
        
        Args:
            data: Input signal data
            window_size: Size of the analysis window
            
        Returns:
            Tuple of (phases, amplitudes) or (None, None) if extraction fails
        """
        if data.shape[0] < window_size:
            print(f'Insufficient data for Hilbert transform: {data.shape[0]} < {window_size}')
            return None, None

        # Extract the analysis window (last window_size samples)
        analysis_window = data[-window_size:]

        # Compute analytic signal using Hilbert transform
        analytic_signal = hilbert(analysis_window)
        phases = np.angle(analytic_signal)
        amplitudes = np.abs(analytic_signal)

        return phases, amplitudes

    def _save_logs(self, filename_prefix: str, path: str = './data') -> None:
        """
        Save timestamp, trigger time, and phase logs to files.
        
        Args:
            filename_prefix: Prefix for the log filenames
            path: Directory path to save the log files
        """
        # Save combined log with all information
        with open(f'{path}/{filename_prefix}_trigger_data.csv', 'w', newline='') as combined_file:
            writer = csv.writer(combined_file)
            writer.writerow(['Timestamp', 'TriggerTime', 'Phase'])
            for ts, tt, phase in zip(self.timestamp_log, self.triggertimes_log, self.phases_log):
                writer.writerow([ts, tt, phase])
        
        # Also save individual files for backward compatibility
        with open(f'{path}/{filename_prefix}_timestamps.csv', 'w', newline='') as ts_file:
            writer = csv.writer(ts_file)
            writer.writerow(['Timestamp'])
            for ts in self.timestamp_log:
                writer.writerow([ts])
        with open(f'{path}/{filename_prefix}_triggertimes.csv', 'w', newline='') as tt_file:
            writer = csv.writer(tt_file)
            writer.writerow(['TriggerTime'])
            for tt in self.triggertimes_log:
                writer.writerow([tt])
        with open(f'{path}/{filename_prefix}_phases.csv', 'w', newline='') as phase_file:
            writer = csv.writer(phase_file)
            writer.writerow(['Phase'])
            for phase in self.phases_log:
                writer.writerow([phase])

        print(f'Logs saved for subject {filename_prefix} in {path}')

test_decider = Decider(num_eeg_channels=64, num_emg_channels=0, sampling_frequency=5000.0)