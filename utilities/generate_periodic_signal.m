function [x, t] = generate_periodic_signal(duration, freq1, freq2)
% GENERATE_PERIODIC_SIGNAL - Generate periodic signal for TSM experiments
%
% Generates the test signal: freq1 + 0.4*freq2
% Default: x = cos(2*pi*261*t) + 0.4*cos(2*pi*783*t)
%
% Input:
%   duration - Signal duration in seconds (default: 1)
%   freq1    - Primary frequency in Hz (default: 261, C4)
%   freq2    - Secondary frequency in Hz (default: 783, G5)
%
% Output:
%   x - Generated signal vector (column vector)
%   t - Time vector (column vector)

% Fixed sampling rate
fs = 44100;

% Set defaults if not provided
if nargin < 1
    duration = 1;
end
if nargin < 2
    freq1 = 261;    % C4 frequency (Hz)
end
if nargin < 3
    freq2 = 783;    % G5 frequency (Hz)
end

% Signal parameters
amp1 = 1.0;     % Primary amplitude
amp2 = 0.4;     % Secondary amplitude

% Generate time vector
t = (0:1/fs:(duration - 1/fs))';

% Generate signal
x = amp1 * cos(2*pi * freq1 * t) + amp2 * cos(2*pi * freq2 * t);

fprintf('Generated periodic signal: %.1f seconds at %d Hz\n', duration, fs);
fprintf('Frequencies: %.0f Hz (amp=%.1f) + %.0f Hz (amp=%.1f)\n', ...
    freq1, amp1, freq2, amp2);

end