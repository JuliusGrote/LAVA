# %%
# packages and paths
import subprocess
import mne
from specparam import SpectralModel
from scipy.signal import firwin
from scipy.io import savemat
from spectrum_analysis import parameterize_spectrum, compute_peak_bands
%matplotlib inline


sub_id = "sub-001"
sess_id = 3
data_load = "../data/sub-001/2025-12-02T162507/"
data_save = "../data/"
rs_file = f"../data/{sub_id}/{sub_id}.set"
coeffs_path =f"../data/{sub_id}_bpfilter.mat"
dec = 10

# %%
matlab_cmd = (
    f"sub_id='{sub_id}'; "
    f"data_load='{data_load}'; "
    f"data_save='{data_save}'; "
    f"sess_id={sess_id}; "
    "run('convert_to_eeglab.m');"
)

result = subprocess.run(
    ["matlab", "-batch", matlab_cmd],
    capture_output=True,
    text=True
)

print(result.stdout)
print(result.stderr)

# %%
# c3 hjorth filter
pre_rs = mne.io.read_raw_eeglab(rs_file, preload=True)

center = "C3"
neighbors = ["FC3", "CP3", "C1", "C5"] 

# get the data
c3 = pre_rs.get_data(picks=center)[::dec, :]
surround = pre_rs.get_data(picks=neighbors)[::dec, :]

c3_hjorth = c3 - surround.mean(axis=0, keepdims=True)

# create a new Raw object for the virtual channel
info = mne.create_info(ch_names=["C3_Hjorth"], sfreq=pre_rs.info["sfreq"]/dec, ch_types=["eeg"])
raw_c3_hjorth = mne.io.RawArray(c3_hjorth, info)
  # Decimate to 500Hz (matches BOSSDevice)
spectrum = raw_c3_hjorth.compute_psd()
psd, freqs = spectrum.get_data(return_freqs=True)

periodic_params, aperiodic_params = parameterize_spectrum([freqs, psd[0,:]])
# compute peak bands
peak = [periodic_params[index][0] for index, param in enumerate(periodic_params) if param[0] < 12 and param[0] > 8][0]
lower, upper = peak - 1, peak + 1
#upper, lower,  = compute_peak_bands(periodic_params, [[8, 12, 0]])

# design FIR filter and retreive coefficients
fs = raw_c3_hjorth.info["sfreq"]
numtaps = 192+1  # Filter order (according to Zrenner et al. 2020) + 1 (must be odd for bandpass)
coefficients = firwin(numtaps, [lower, upper], fs=fs, pass_zero='bandpass')

# save coefficients to .mat file
savemat(coeffs_path.format(sub=sub_id), {'coefficients': coefficients})


