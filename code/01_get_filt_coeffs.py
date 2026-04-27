# %%
# packages and paths
import subprocess
from pathlib import Path
import mne
from scipy.signal import firwin
from scipy.io import savemat
from spectrum_analysis import parameterize_spectrum

sub_id = "sub-001"
sess_id = 3
data_load = fr"C:\ProgramData\Bittium Biosignals Ltd\NeurOne64\SessionData\NeurOneUser\LAVA_{sub_id[-3:]}\2025-12-02T162507"
print(f"data_load: {data_load}")
script_dir = Path(__file__).resolve().parent
repo_root = script_dir.parent
data_save = repo_root / "data"
subject_dir = data_save / sub_id
rs_file = subject_dir / f"{sub_id}.set"
coeffs_path = data_save / f"{sub_id}_bpfilter.mat"
dec = 10
print(f"root: {repo_root}")
subject_dir.mkdir(parents=True, exist_ok=True)

# %%
matlab_cmd = (
    f"root='{repo_root}'; "
    f"sub_id='{sub_id}'; "
    f"data_load='{data_load}'; "
    f"data_save='{data_save.as_posix()}'; "
    f"sess_id={sess_id}; "
    "run('convert_to_eeglab.m');"
)

result = subprocess.run(
    ["matlab", "-batch", matlab_cmd],
    capture_output=True,
    text=True,
    cwd=str(script_dir)
)

print(result.stdout)
print(result.stderr)

if result.returncode != 0:
    raise RuntimeError(
        f"MATLAB conversion failed with exit code {result.returncode}.\n"
        f"stdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}"
    )

if not rs_file.exists():
    raise FileNotFoundError(
        f"Expected EEGLAB file was not created: {rs_file}.\n"
        "Check MATLAB conversion logs above and confirm input data path/session are correct."
    )

# %%
# c3 hjorth filter
pre_rs = mne.io.read_raw_eeglab(str(rs_file), preload=True)

center = "C3"
neighbors = ["FC3", "CP3", "C1", "C5"] 

# get the data
c3 = pre_rs.get_data(picks=center)[:, ::dec]
surround = pre_rs.get_data(picks=neighbors)[:, ::dec]

c3_hjorth = c3 - surround.mean(axis=0, keepdims=True)

# create a new Raw object for the virtual channel
info = mne.create_info(ch_names=["C3_Hjorth"], sfreq=pre_rs.info["sfreq"]/dec, ch_types=["eeg"])
raw_c3_hjorth = mne.io.RawArray(c3_hjorth, info)
  # Decimate to 500Hz (matches BOSSDevice)
spectrum = raw_c3_hjorth.compute_psd()
psd, freqs = spectrum.get_data(return_freqs=True)

periodic_params, aperiodic_params = parameterize_spectrum([freqs, psd[0,:]], save_fig=repo_root / "figures" / f"{sub_id}_psd.png")
# compute peak bands
peak = [periodic_params[index][0] for index, param in enumerate(periodic_params) if param[0] < 12 and param[0] > 8][0]
lower, upper = peak - 1, peak + 1
#upper, lower,  = compute_peak_bands(periodic_params, [[8, 12, 0]])

# design FIR filter and retrieve coefficients
fs = raw_c3_hjorth.info["sfreq"]
numtaps = 80+1  # Filter order (according to Zrenner et al. 2020) + 1 (must be odd for bandpass)

print(f"Designing FIR filter with fs={fs}, numtaps={numtaps}, lower={lower}, upper={upper}")

coefficients = firwin(numtaps, [lower, upper], fs=fs, pass_zero='bandpass')

# save coefficients to .mat file
savemat(str(coeffs_path), {'coefficients': coefficients})


