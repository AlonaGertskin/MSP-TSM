function config = project_config(experiment)
% PROJECT_CONFIG - Configuration parameters for TSM project
%
% Usage:
%   config = project_config()           % Returns global config
%   config = project_config('exp1')     % Returns experiment 1 config

if nargin == 0
    % Global parameters shared across all experiments
    config.fs = 44100;                % Default sampling rate (Hz)
    config.output_dir = 'outputs';    % Output directory
else
    switch experiment
        case 'exp1'
            config = config_exp1();
        case 'exp1_sweep'
            config = config_exp1_sweep();
        otherwise
            error('Unknown experiment: %s', experiment);
    end
end

end

function config = config_exp1()
% Configuration for Experiment 1: Pitch Modification with Periodic Signals

% Get global config first
config = project_config();

% Add experiment-specific parameters
config.duration = 1;                  % Signal duration (seconds)
config.freq1 = 261;                   % C4 frequency (Hz)
config.freq2 = 783;                   % G5 frequency (Hz) - 3rd harmonic
config.amp1 = 1.0;                    % Primary amplitude
config.amp2 = 0.4;                    % Secondary amplitude
config.target_g4 = 392;               % G4 (Dominant)
config.target_f3 = 175;               % F3 (Sub-Dominant)

end

function config = config_exp1_sweep()
% Configuration for Experiment 1 parameter sweep

% Parameter sweep values for OLA
% config.ola_synHop_values = [512,1024, 2048];
% config.ola_winLen_values = [256, 512];
% config.ola_winType_values = [2, 4];      % 1=sine, 2=Hann
config.ola_synHop_values = [1024];
config.ola_winLen_values = [256];
config.ola_winType_values = [2];      % 1=sine, 2=Hann
%Parameter sweep values for WSOLA
% config.wsola_tolerance_values = [512, 1024, 2048, 4096];
config.wsola_tolerance_values = [512];

end