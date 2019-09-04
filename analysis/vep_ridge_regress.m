% ------------------------------------------------------------------------
%                 Ridge Regression for Decoding cVEP
% -------------------------------------------------------------------------
function vep_ridge_regress(epochedEEG,stim_seq)

STIM = 1;
CONT = 3;
CH = [12,19]; % [5,10:14,18:20]; % [12,18:20]; % [5,9:15,17:21,24:26]; %[12,17:21,24:26];
REFRESH = 60;
SRATE = 512;

% define time window of EEG response to each code or edge
segWindow = round([-150, 50] * SRATE / 1000);    % samples

% extract timing of code sequence or that of rising / falling edges
on_fmc = cell(1,2);
off_fmc = cell(1,2);
rise_fmc = cell(1,2);
fall_fmc = cell(1,2);
flat_fmc = cell(1,2);
for loc_it = 1:2
    on_fmc{loc_it}  = find(stim_seq.fmc{loc_it} == 1);
    off_fmc{loc_it} = find(stim_seq.fmc{loc_it} == 0);
    rise_fmc{loc_it} = find((stim_seq.fmc{loc_it} - [0 stim_seq.fmc{loc_it}(1:end-1)]) == 1);
    fall_fmc{loc_it} = find((stim_seq.fmc{loc_it} - [0 stim_seq.fmc{loc_it}(1:end-1)]) == -1);
    flat_fmc{loc_it} = find((stim_seq.fmc{loc_it} - [0 stim_seq.fmc{loc_it}(1:end-1)]) == 0);
end

% extract EEG epochs time-locked to corresponding code sequence
X_on = []; Y_on = []; X_off = []; Y_off = [];
X_rise = []; Y_rise = []; X_fall = []; Y_fall = []; X_flat = []; Y_flat = [];
for trial_id = 1:10
    for LOC = 1:2
        data = epochedEEG{STIM,CONT,LOC,trial_id}(CH,:);
        
        [X,Y] = epoch_eeg_code(data,on_fmc{LOC},1,segWindow,SRATE,REFRESH);
        X_on = [X_on; X]; Y_on = [Y_on; Y];
        [X,Y] = epoch_eeg_code(data,off_fmc{LOC},-1,segWindow,SRATE,REFRESH);
        X_off = [X_off; X]; Y_off = [Y_off; Y];
        [X,Y] = epoch_eeg_code(data,rise_fmc{LOC},1,segWindow,SRATE,REFRESH);
        X_rise = [X_rise; X]; Y_rise = [Y_rise; Y];
        [X,Y] = epoch_eeg_code(data,fall_fmc{LOC},-1,segWindow,SRATE,REFRESH);
        X_fall = [X_fall; X]; Y_fall = [Y_fall; Y];
        
    end
end

% separate training and testing data
CV_FOLD = 10;
alpha_list = 10.^(0:8);

% select conditions
% X1 = X_on; X0 = X_off;
% Y1 = Y_on; Y0 = Y_off;
X1 = X_rise; X0 = X_fall;
Y1 = Y_rise; Y0 = Y_fall;

% obtain number of trials for both codes {1,0}
[total_trial_1, feat_size] = size(X1);
[total_trial_0, ~] = size(X0);
block_size_1 = floor(total_trial_1 / CV_FOLD);
block_size_0 = floor(total_trial_0 / CV_FOLD);

acc_ls = zeros(1,CV_FOLD);
acc_ridge = zeros(length(alpha_list),CV_FOLD);

% cross validation 
for cv_id = 1:CV_FOLD
    
    % divide test and training trials for each condition
    test_trial_1_id = (cv_id-1)*block_size_1+1 : cv_id*block_size_1;
    test_trial_0_id = (cv_id-1)*block_size_0+1 : cv_id*block_size_0;
    train_trial_1_id = setdiff(1:total_trial_1, test_trial_1_id);
    train_trial_0_id = setdiff(1:total_trial_0, test_trial_0_id);

    X_test = [X1(test_trial_1_id,:); X0(test_trial_0_id,:)];
    Y_test = [Y1(test_trial_1_id,:); Y0(test_trial_0_id,:)];
    X_train = [X1(train_trial_1_id,:); X0(train_trial_0_id,:)];
    Y_train = [Y1(train_trial_1_id,:); Y0(train_trial_0_id,:)];
    
    % linear regression and least-square solution
    b_ls = (X_train' * X_train) \ (X_train' * Y_train);
    Y_test_est_ls = 2 * (X_test * b_ls > 0) - 1;
    acc_ls(cv_id) = sum(Y_test_est_ls == Y_test) / length(Y_test);
    
    % ridge regression
    for alpha_id = 1:length(alpha_list)
        alpha = alpha_list(alpha_id);
        b_ridge = (X_train' * X_train + alpha * eye(feat_size) ) \ (X_train' * Y_train);
        Y_test_est_ridge = 2 * (X_test * b_ridge > 0) - 1;
        acc_ridge(alpha_id,cv_id) = sum(Y_test_est_ridge == Y_test) / length(Y_test);
    end
    
end

disp('Linear Regression Accuracy:')
disp(mean(acc_ls))
disp('Ridge Regression Accuracy:')
disp(max(mean(acc_ridge,2)))

end



function [X,Y] = epoch_eeg_code(data,code,label,segWindow,SRATE,REFRESH)

% extract EEG epochs time-locked to code label
X = []; Y = [];
for it = 1:length(code)
    startIndex = round((code(it)-1) / REFRESH * SRATE) - segWindow(1);
    if startIndex+segWindow(2) <= length(data)
        temp_epoch = data(:, startIndex+1+segWindow(1) : startIndex+segWindow(2));
        X = [X; temp_epoch(:)' ];
        Y = [Y; label];
    end
end

end