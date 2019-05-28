Data Description:

Note: Run vep_analysis_script() as an example for loading and analyzing data using deconvolution, ridge regression, and TRCA. 


A. Epoched EEG Data:

1. First Subject: S1_epochedEEG_exp1.mat
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


2. Second and Third Subjects: 

a. Experiment 1: S2_epochedEEG_exp1.mat & S3_epochedEEG_exp1.mat
% Variable name: epochedEEG_exp1
%   - 4-dim cell structure: 5 x 3 x 2 x 10
%           5 (stimulation patterns: 1: fmc: 30-Hz m-sequence, 2: standard m-sequence, 3: 30-Hz SSVEP, 4: fmc on color image, 5: no flickering)
%           3 (contrast levels: 2, 4, 8)
%           2 (locations: 1: left, 2: right half of the screen)
%           10 (trial ID, total 10 trials for each condition)
%   - Each cell contains epoched EEG data 33 x 1742
%           33 (S2) / 31 (S3)  (number of channels, excluding 1 / 3 bad channel)
%           1742 (epoched data samples from -0.2 to 3.2 sec time locked to
%           the onset of stimulation, with sampling rate 512)

b. Experiment 2: S2_epochedEEG_exp2.mat & S3_epochedEEG_exp2.mat
% Variable name: epochedEEG_exp2
%   - 2-dim cell structure: 11 x 10
%           11 (locations: 20, 30, 35, 40, 45, 50, 55, 60, 65, 70, 80 % of the screen, <50% - left screen, >50% - right screen)
%           10 (trial ID, total 10 trials for each condition)
%   - Each cell contains epoched EEG data 33 x 1742
%           33 (S2) / 31 (S3) (number of channels, excluding 1 / 3 bad channel)
%           1742 (epoched data samples from -0.2 to 3.2 sec time locked to
%           the onset of stimulation, with sampling rate 512)



B. Stimulation Sequence Data: S1_stim_data.mat or S1_stim_data.mat
% Variable name: stim_data
%   - 4-dim cell structure: 4 x 3 x 2 x 10 (same as epochedEEG)
%   - Each cell contains stimuli sequence (2 x 1640) that cooresponds to the epochedEEG trial
%           2 (1: raw sequence {1,0,-1} interpolated to 512 Hz with zeros in between 60Hz frames
%             (2: edge sequence {1,0,-1} rising edge and falling edge
%           1640 time points

% Variable name: stim_seq
%   - struct with three fields
%       .fmc (30Hz frequency-modulated msequence)
%       .fmc_inv (fmc with inverted code for left and right screen)
%       .mseq (msequence)
%   - each field consists of 1x2 cell array
%       {1}: 1x180 (60Hz refresh rate x 3 sec) code sequence of left screen
%       {2}: 1x180 (60Hz refresh rate x 3 sec) code sequence of right screen
%   - Note: there is a constant delay between EEG and code-sequence streams


