% ------------------------------------------------------------------------
%               TRCA-based VEP detection algorithm
% -------------------------------------------------------------------------
function [mean_acc, mu_ci] = vep_trca(epochedEEG,STIM,CONT)

fprintf('Results of the ensemble TRCA-based method.\n');

% select stimulation pattern:
numTrial = 10;
SRATE = 512;
plot_fig = 0;

% define parameters
CH = 1:size(epochedEEG{1,1,1,1},1); %[5,10:14,18:20]; % {[12,19], [12,18:20], [5,10:14,18:20], [5,9:15,17:21,24:26], [1:15,17:21,24:26], [1:32]};
data_len_list = 3.4;   % sec
data_offset = 0;   % sec
is_ensemble = 1;    % use ensemble classifier (one classifier for each code seq)

% define conditions for testing
var_list = data_len_list;    % different channel combinations

mean_acc = zeros(1, length(var_list));
mu_ci = zeros(length(var_list),2);

% main loop - cross validation for each condition
for var_i = 1:length(var_list)
    
    data_len = var_list(var_i);
    data_range = floor(data_offset*SRATE) + (1:1:floor(data_len*SRATE));
    NSAMP = length(data_range);
    
    % extract EEG epochs
    trial_eeg = zeros(2,length(CH),NSAMP,numTrial);
    for loc_i = 1:2
        for tr_i = 1:10
            trial_eeg(loc_i,:,:,tr_i) = epochedEEG{STIM,CONT,loc_i,tr_i}(CH,data_range);
        end
    end
    
    % Leave-one-trial-out cross validation classification accuracy
    labels = [1, 2];
    accs = zeros(1,numTrial);
    for cv_i = 1:1:numTrial
        
        % Training stage
        traindata = trial_eeg;
        traindata(:, :, :, cv_i) = [];
        model = train_trca(traindata);
        
        % Test stage
        testdata = squeeze(trial_eeg(:, :, :, cv_i));
        estimated = test_trca(testdata, model, is_ensemble);
        
        % Evaluation
        is_correct = (estimated==labels);
        accs(cv_i) = mean(is_correct)*100;
        
    end % loocv_i
    
    % Summarize results
    alpha_ci = 0.05;    % alpha for computing confidence interval
    ci = 1-alpha_ci;    % confidence interval
    [mu, ~, muci, ~] = normfit(accs, alpha_ci);
    fprintf('Mean accuracy = %2.2f %% (%2d%% CI: %2.2f - %2.2f %%)\n',...
        mu, ci, muci(1), muci(2));
    
    mean_acc(var_i) = mu;
    mu_ci(var_i,:) = muci;
    
end

% plot classification results vs. number of channels for each contrast
if plot_fig
    figure, plot(var_list, mean_acc,'linewidth',2);
    xlabel('Training data legnth (sec)'); ylabel('Cross validation accuracy (''%)');
    set(gca,'fontsize',12)
    ylim([0 100]);
end

end
