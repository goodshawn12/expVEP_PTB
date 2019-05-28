addpath('data');

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
tmp = load('S3_epochedEEG_exp1.mat');
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
