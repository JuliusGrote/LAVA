function coeffs = create_filter_coeffs()
    % Create bandpass filter coefficients matching LAVA_NeuroSimo.py firwin parameters
    %
    % Python equivalent:
    % firwin(filter_order + 1, [low_normalized, high_normalized], 
    %        pass_zero=False, window='hamming')
    %
    % Parameters from Python:
    % - Filter order: 128 (129 coefficients)
    % - Low cutoff: 8 Hz
    % - High cutoff: 14 Hz
    % - Sampling rate: 1000 Hz
    % - Window: Hamming
    
    % Filter parameters
    filter_order = 128;
    low_cutoff = 8.0;    % Hz
    high_cutoff = 14.0;  % Hz
    sampling_rate = 5000.0;  % Hz
    
    % Normalize cutoff frequencies to [0, 1] where 1 is Nyquist frequency
    nyquist = sampling_rate / 2;
    normalized_cutoffs = [low_cutoff, high_cutoff] / nyquist;
    
    % Design FIR bandpass filter using designfilt
    % This creates the same filter as Python's firwin with pass_zero=False
    d = firls(80, [0 (11 + [-5 -2 +2 +5]) (500/2)]/(500/2), [0 0 1 1 0 0], [1 1 1] );
    
    % firls returns coefficients directly, not a filter object
    coeffs = d;
    
    % Display filter information
    fprintf('Generated bandpass filter coefficients with order %d for %.1f-%.1f Hz band.\n', ...
        filter_order, low_cutoff, high_cutoff);
    fprintf('Number of coefficients: %d\n', length(coeffs));
    fprintf('Filter Coefficients:\n');
    
    % Display all digits using format long
    format long
    disp(coeffs');
    
    % Save coefficients to .mat file for direct loading in Python
    save('filter_coeffs.mat', 'coeffs');

    % Optional: visualize frequency response
    % fvtool(d);
end
