function vep_analysis_script()

% load epoched EEG responses to visual stimulation

% Variable name: epochedEEG
%   - 4-dim cell structure: 4 x 3 x 2 x 10
%           4 (stimulation patterns: 1: fmc-VEP (independent) 2: fmc-VEP (inverted) 3: fmc-VEP (chrome) 4: m-sequence)
%           3 (contrast levels: 4, 8, 64)
%           2 (locations: 1: left, 2: right half of the screen)
%           10 (trial ID, total 10 trials for each condition)
%   - Each cell contains epoched EEG data 34 x 1640
%           34 (number of channels)
%           1640 (epoched data samples from -0.2 to 3 sec time locked to
%           the onset of stimulation, with sampling rate 512)
load('S1_epochedEEG_exp1.mat');

% Variable name: stim_data
%   - 4-dim cell structure: 4 x 3 x 2 x 10 (same as epochedEEG)
%   - Each cell contains stimuli sequence (2 x 1640) that cooresponds to the epochedEEG trial
%           2 (1: raw sequence {1,0,-1} interpolated to 512 Hz with zeros in between 60Hz frames
%             (2: edge sequence {1,0,-1} rising edge and falling edge
%           1640 time points
load('S1_stim_data.mat');

% Variable name: stim_seq
%   - struct with three fields
%       .fmc (30Hz frequency-modulated msequence)
%       .fmc_inv (fmc with inverted code for left and right screen)
%       .mseq (msequence)
%   - each field consists of 1x2 cell array
%       {1}: 1x180 (60Hz refresh rate x 3 sec) code sequence of left screen
%       {2}: 1x180 (60Hz refresh rate x 3 sec) code sequence of right screen
%   - Note: there is a constant delay between EEG and code-sequence streams
load('S1_stim_seq.mat');

% deconvolution analysis
vep_deconv(epochedEEG,stim_data)

% % ridge regression analysis
vep_ridge_regress(epochedEEG,stim_seq)

% task-related component analysis
vep_trca(epochedEEG)


end





