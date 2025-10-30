function out = createTimeAxis(timeWindow,samplingFrequency)

out = (timeWindow(1)+1/samplingFrequency*1000):1/samplingFrequency*1000:timeWindow(2); % time axis (ms)