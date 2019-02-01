% Clear the workspace
close all;
clearvars;
sca;

% add LSL library into path (for marker stream output)
addpath(genpath('liblsl-Matlab'))


%% ------------------------------------------------------------------------
%              Define Configurations of the Experiment
%--------------------------------------------------------------------------
config = [];

% Stimuli modes - code [A B]
%   A: Stimulation code     - 0: frequency-modulated code 1: m-sequence 2: ssvep, 3: image, 4: text
%   B: Modulation modes     - [A==0] 0: normal, 1: opposite codes, 2: chromatical modulation, 3: slow (15Hz) modulation 

% Examples (all available options):
%   'fmc_normal':           [0 0] 
%   'fmc_opposite':         [0 1]
%   'fmc_chrome':           [0 2]
%   'fmc_slow':             [0 3]
%   'mseq':                 [1 0]
%   'ssvep':                [2 0]
%   'fmc_single_image':     [3 0]   (for demo only)
%   'fmc_single_text':      [4 0]   (for demo only)

% setting - experiment tasl 1
config.RUN_EXP_1 = 1;        % run the experiment task 1 - contrast levels and stimuli modes
config.MODE = {[0 0], [0 1], [1 0], [0 2]};     % select all the stimuli modes for the experiment
config.CONTRAST_LIST = [4, 8, 64];  % conditions for experiment task 1
config.LIST_EACH_MODE = {1:3, 1:3, 1:3, 1:2};   % contrast list index for each mode
config.NUM_TRIAL_CONTRAST = 10;     % number of trials per contrast level per location (left/right)
config.BREAK_INTERVAL = 22;        % break every N trials

% setting - experiment tasl 2
config.RUN_EXP_2 = 0;        % run the experiment task 2 - fixation location
config.MODE_LOC = [0,0];    % experiment 2
config.NUM_TRIAL_LOCATION = 2;     % number of trials per contrast level per location (left/right)
config.LOCATION_LIST = [25, 40, 60, 75];    % conditions for experiment task 2
config.CONTRAST = 30;        % [0 127]

% setting - demo experiment stimuli
config.RUN_DEMO = 0;         % run demo - show stimuli only
config.MODE_DEMO = [0,2];   % demo 
config.NUM_REPEAT_DEMO = 5;

% other experiment settings
config.DEMO = 0;            % skip synchronization test if in demo mode
config.ENABLE_LSL = 1;      % enable sending LSL marker stream 


% screen setting
config.FULLSCREEN = 1;      % if full screen fails, try specify the window size below
config.WIN_WIDTH = 1600;    % specify the screen resolution if not full screen (in pixel)
config.WIN_HEIGHT = 1200;   
config.REFRESH = 60;     % refresh rate in Hz

% contrast level between 
config.STIM_LEN = 3;        % second
config.BASELINE = 128;      % [0 255]

% smooth boundary
config.SMOOTH = 1;          % 0: no smooth, 1: linear smooth, 2: probabilistic smooth
config.SMOOTH_WIDTH = 100;  % smoothing width in pixel

% specify image and text files
config.filename_image = 'img/img01.jpg';
config.filename_text = 'img/text00.png';


%% run the experiment and record response matrix
respMat = main(config);
save('respMat.mat','respMat')
