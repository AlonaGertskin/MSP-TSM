function y = ola(x, alpha, synHop, winLen)
% OLA - Overlap-Add time-scale modification
%
% Simple wrapper for wsolaTSM.m with tolerance = 0
%
% Input:
%   x      - Input signal
%   alpha  - Stretch factor (Hs/Ha)
%   synHop - Synthesis hop size (optional, default: 128)
%   winLen - Window length (optional, default: 256)
%
% Output:
%   y - Time-scale modified signal

% Set default parameters
if nargin < 3 || isempty(synHop)
    synHop = 128;
end
if nargin < 4 || isempty(winLen)
    winLen = 256;
end

% Create parameter structure (following the example)
paramOLA.tolerance = 0;
paramOLA.synHop = synHop;
paramOLA.win = win(winLen, 2);  % Using the TSM toolbox win function

% Call wsolaTSM with parameter structure
y = wsolaTSM(x, alpha, paramOLA);

end