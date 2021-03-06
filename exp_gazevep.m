
% Clear the workspace
close all;
clearvars;
sca;

% add LSL library into path (for marker stream output)
addpath(genpath('C:\Users\shawn\Documents\exp_vep\liblsl-Matlab'))
addpath('img')

%% ------------------------------------------------------------------------
%              Define Configurations of. the Experiment
%--------------------------------------------------------------------------
config = [];

% screen setting
config.FULLSCREEN = 0;      % if full screen fails, try specify the window size below
config.WIN_WIDTH = 1280;    % specify the screen resolution if not full screen (in pixel)
config.WIN_HEIGHT = 960;   
config.REFRESH = 60;     % refresh rate in Hz

% contrast level between 
config.STIM_LEN = 3;        % second
config.BASELINE = 128;      % [0 255]

% smooth boundary
config.SMOOTH = 1;          % 0: no smooth, 1: linear smooth, 2: probabilistic smooth
config.SMOOTH_WIDTH = 100;  % smoothing width in pixel

% other experiment settings
config.DEMO = 0;            % skip synchronization test if in demo mode
config.ENABLE_LSL = 1;      % enable sending LSL marker stream 

% Stimuli modes:
%   1: 'fmc' (gray)
%   2: 'mseq' (gray)
%   3: 'ssvep' (gray)
%   4: 'fmc_image' (image)

% setting - experiment task 1
config.RUN_EXP_1 = 1;        % run the experiment task 1 - contrast levels and stimuli modes
config.MODE = [1, 2, 3, 4];     % select all the stimuli modes for the experiment 
config.CONTRAST_LIST = [2, 8, 16];  % conditions for experiment task 1
config.NUM_TRIAL_CONTRAST = 10;     % number of trials per contrast level per location (left/right)
config.BREAK_INTERVAL = 12;        % break every N trials
config.SESS_INTERNAL = 5;         % break every M blocks
config.NUM_IMG = 10;

% setting - calibration - blinking task for event synchronization
config.RUN_CALIB = 1;

% setting - demo experiment stimuli
config.RUN_DEMO = 0;         % run demo - show stimuli only
config.MODE_DEMO = 4;   % demo 
config.DEMO_CONTRAST_LEVEL = 2;
config.DEMO_IMAGE = 1;
config.DEMO_LOCATION = 1;
config.NUM_REPEAT_DEMO = 5;


%% run the experiment and record response matrix
respMat = main(config);
save(sprintf('respMat_%d.mat',floor(rand*10000)),'respMat')
