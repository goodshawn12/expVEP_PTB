% Clear the workspace
close all;
clearvars;
sca;

% add LSL library into path (for marker stream output)
% addpath(genpath('C:\Users\shawn\Desktop\NSF_BIC\liblsl-Matlab'))


%% ------------------------------------------------------------------------
%              Define Configurations of the Experiment
%--------------------------------------------------------------------------
config = [];

% Stimulation modes: 'A_B_C'
%   A: Stimulation code  {'fmc': frequency-modulated code, 
%                         'mseq': m-sequence 
%                         'ssvep': steady-state visual-evoked potential}
%   B: Split screens     {'binary': two independent codes,
%                         'opposite': opposite codes,
%                         'single': single stimuli, for image/text only}
%   C: Miscellaneous     {'': luminance modulation, 
%                         'chrome': chromatic modulation,
%                         'slow': 15 Hz modulation,
%                         'image': image (single screen only).
%                         'text': text (single screen only)}
% Examples (all available options)
%   'fmc_binary', 'fmc_binary_chrome', 'fmc_binary_slow'
%   'fmc_opposite', 'fmc_opposite_chrome', 'fmc_opposite_slow'
%   'mseq_opposite, 'mseq_binary, 
%   'ssvep'
%   'fmc_single_image', 'fmc_single_text'
config.MODE = 'ssvep'; 


% experiment setting
config.RUN_EXP_1 = 0;        % run the experiment task 1 - contrast levels
config.RUN_EXP_2 = 0;        % run the experiment task 2 - fixation location
config.RUN_DEMO = 1;         % run demo - show stimuli only

%
config.NUM_TRIAL_CONTRAST = 1;     % number of trials per contrast level per location (left/right)
config.NUM_TRIAL_LOCATION = 2;     % number of trials per contrast level per location (left/right)
config.STIM_LEN = 2;        % second
config.DEMO = 1;            % skip synchronization test if in demo mode
config.ENABLE_LSL = 0;      % enable sending LSL marker stream 

% screen setting
config.FULLSCREEN = 1;      % if full screen fails, try specify the window size below
config.WIN_WIDTH = 1920;    % specify the screen resolution if not full screen (in pixel)
config.WIN_HEIGHT = 1080;   
config.REFRESH = 60;     % refresh rate in Hz

% contrast level between 
config.BASELINE = 128;
config.CONTRAST = 5;     
config.CONTRAST_LIST = [1, 4, 16, 64, 96];  % conditions for experiment task 1
config.LOCATION_LIST = [25, 40, 60, 75];    % conditions for experiment task 2

% smooth boundary
config.SMOOTH = 1;          % smooth at boundary of left and right screens
config.SMOOTH_WIDTH = 100;  % smoothing width in pixel

% specify image and text files
config.filename_image = 'img/img01.jpg';
config.filename_text = 'img/text00.png';


%% run the experiment and record response matrix
respMat = main(config);

