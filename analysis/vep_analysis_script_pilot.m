addpath('data');

% a. Experiment 1: S2_epochedEEG_exp1.mat & S3_epochedEEG_exp1.mat
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
%
% b. Experiment 2: S2_epochedEEG_exp2.mat & S3_epochedEEG_exp2.mat
% Variable name: epochedEEG_exp2
%   - 2-dim cell structure: 11 x 10
%           11 (locations: 20, 30, 35, 40, 45, 50, 55, 60, 65, 70, 80 % of the screen, <50% - left screen, >50% - right screen)
%           10 (trial ID, total 10 trials for each condition)
%   - Each cell contains epoched EEG data 33 x 1742
%           33 (S2) / 31 (S3) (number of channels, excluding 1 / 3 bad channel)
%           1742 (epoched data samples from -0.2 to 3.2 sec time locked to
%           the onset of stimulation, with sampling rate 512)

tmp = load('S2_epochedEEG_exp1.mat');
epochedEEG = tmp.epochedEEG_exp1;

% task-related component analysis
mean_acc = zeros(4,3);
for STIM = 1:4
    for CONT = 1:3
        [mean_acc(STIM,CONT), ~] = vep_trca(epochedEEG,STIM,CONT);
    end
end

figure, plot(mean_acc','linewidth',2);
legend('FMC','MSEQ','SSVEP','FMC-IMG');
xlabel('Contrast'); set(gca,'XTick',[1,2,3],'XTickLabel',{'2','4','8'},'fontsize',14);
ylabel('Cross Validation Accuracy'); ylim([0 100]);


% load behavioral results
tmp_beh = load('respMat_051519S3.mat');
respMat = tmp_beh.respMat{2};
vep_behavior(respMat);

%{
% deconvolution analysis
vep_deconv(epochedEEG,stim_data)

% % ridge regression analysis
vep_ridge_regress(epochedEEG,stim_seq)
%}
