% Clear the workspace
close all;
clearvars;
sca;

% add LSL library into path
addpath(genpath('C:\Users\shawn\Desktop\NSF_BIC\liblsl-Matlab'))

%--------------------------------------------------------------------------
%              Define Configurations of the Experiment
%--------------------------------------------------------------------------

config = [];


% stimuli modes
%   'fmc_opposite', 'fmc_opposite_chrome', 'fmc_opposite_slow'
%   'fmc_binary', 'fmc_binary_chrome', 'fmc_binary_slow'
%   'mseq_opposite, 'mseq_binary, 'ssvep'
%   'fmc_single_image', 'fmc_single_text'
config.MODE = 'fmc_binary'; 
config.MODE_STIM = 0;   % 0: frequency-modulated code 1: m-sequence 2: ssvep
config.MODE_SIDE = 0;   % 0: two independent codes, 1: opposite codes, (3: single stimuli, for image/text only)
config.MODE_COLOR = 0;  % 0: luminance modulation, 1: chromatic modulation
config.MODE_NORM = 1;   % 0: slow (15 Hz) modulation, 1: normal (30 Hz) modulation

% experiment setting
config.RUN_EXP_1 = 1;        % run the experiment task 1 - contrast levels
config.RUN_EXP_2 = 1;        % run the experiment task 2 - fixation location
config.RUN_DEMO = 0;         % run demo - show stimuli only

%
config.NUM_TRIAL_CONTRAST = 1;     % number of trials per contrast level per location (left/right)
config.NUM_TRIAL_LOCATION = 2;     % number of trials per contrast level per location (left/right)
config.STIM_LEN = 2;        % second
config.DEMO = 1;            % skip synchronization test if in demo mode
config.ENABLE_LSL = 1;      % enable sending LSL marker stream 

% screen setting
config.FULLSCREEN = 1;      % if full screen fails, try specify the window size below
config.WIN_WIDTH = 1920;    % specify the screen resolution if not full screen (in pixel)
config.WIN_HEIGHT = 1080;   
config.REFRESH = 60;     % refresh rate in Hz

% contrast level between 
config.BASELINE = 128;
config.CONTRAST = 5;     
config.CONTRAST_LIST = [1, 4, 16, 64, 96];
config.LOCATION_LIST = [25, 40, 60, 75];

% smooth boundary
config.SMOOTH = 1;
config.SMOOTH_WIDTH = 100;  % pixel

% specify image and text files
config.filename_image = 'img/img01.jpg';
config.filename_text = 'img/text00.png';

respMat = main(config);


%% TODO
%{
1. Stimuli presentation
a. SSVEP binary (done)
b. m-sequence binary (done)
    - automatic select num_bit by defined code_length
    - fix code: using the same initial states (done)
c. chromatic modulation (done)
d. 15Hz modulation (done)
e. smooth boundary (done)

%}