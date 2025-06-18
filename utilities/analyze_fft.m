function [freq_axis, magnitude, peaks] = analyze_fft(x, fs, plot_result, freq_range)

% ANALYZE_FFT - Perform FFT analysis on signal
%
% Input:
%   x           - Input signal (column vector)
%   fs          - Sampling frequency (default: 44100)
%   plot_result - Boolean to plot result (default: false)
%
% Output:
%   freq_axis   - Frequency axis (Hz)
%   magnitude   - Magnitude spectrum
%   peaks       - Structure with peak frequencies and magnitudes

% Set defaults
if nargin < 2
    fs = 44100;
end
if nargin < 3
    plot_result = false;
end

if nargin < 4
    freq_range = [0, 2000];  % Default range for backward compatibility
end

% Ensure signal is column vector
if size(x, 1) == 1
    x = x(:);
end

% FFT parameters
N = length(x);
NFFT = 2^nextpow2(N);  % Zero-pad to next power of 2 for efficiency

% Compute FFT
X = fft(x, NFFT);
magnitude = abs(X);

% Create frequency axis (only positive frequencies)
freq_axis = (0:NFFT/2-1) * fs / NFFT;
magnitude = magnitude(1:NFFT/2);

% Find peaks (look for frequencies up to 2kHz)
max_freq_idx = find(freq_axis <= 2000, 1, 'last');
search_mag = magnitude(1:max_freq_idx);
search_freq = freq_axis(1:max_freq_idx);

% Improved peak detection
min_peak_height = max(search_mag) * 0.2;  % Peaks must be at least 20% of max
min_peak_distance = round(length(search_mag) * 0.01);  % Minimum distance between peaks

[peak_mags, peak_locs] = findpeaks(search_mag, ...
    'MinPeakHeight', min_peak_height, ...
    'MinPeakDistance', min_peak_distance, ...
    'SortStr', 'descend');

% Get peak frequencies
peak_freqs = search_freq(peak_locs);

% Keep only top 3 most significant peaks
num_peaks = min(3, length(peak_freqs));
peak_freqs = peak_freqs(1:num_peaks);
peak_mags = peak_mags(1:num_peaks);

% Store results
peaks.frequencies = peak_freqs;
peaks.magnitudes = peak_mags;

% Print peaks
fprintf('Top frequency peaks:\n');
for i = 1:num_peaks
    fprintf('  %.1f Hz (magnitude: %.2f)\n', peak_freqs(i), peak_mags(i));
end

% Plot if requested
if plot_result
    fig = figure('Visible', 'on');  % Show plot
else
    fig = figure('Visible', 'off'); % Create plot but don't show
end

plot(freq_axis, magnitude, 'b-', 'LineWidth', 1);
hold on;
if ~isempty(peak_freqs)
    plot(peak_freqs, peak_mags, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
end
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('FFT Magnitude Spectrum');
xlim(freq_range);  % Show up to 2kHz for audio signals
legend('Spectrum', 'Peaks', 'Location', 'best');

end