
% addpath EEGLAB
cd('C:\Users\shawn\Desktop\NSF_BIC\expVEP_PTB')


%% ------------------------------------------------------------------------
%               Load EEG dataset and extract EEG epochs
% -------------------------------------------------------------------------

% load preprocessed EEG data
filename = 'expVEP_013119_EEG.set';
filepath = 'C:\Users\shawn\Desktop\NSF_BIC\expVEP_PTB\data\';
EEG = pop_loadset([filepath filename]);

% load behavioral results and stimulation codes
load('data\respMat_013119_pilot.mat')
code = respMat{1};

% extract EEG epoch
epochWindow = [-0.2, 3];   % define time window for EEG epoch (in sec)
numTrial = 10;
expression = 'Stim(?<stim>\d+)_CT(?<cont>\d+)_Loc(?<loc>\d+)';  % decode event markers
stimList = {'00','01','02','10'};   % stimulation pattern: 00: fmc-VEP (independent) 01: fmc-VEP (inverted) 02: fmc-VEP (chrome) 10: mseq
contList = {'4','8','64'};          % contrast level: +- 4, 8, 64
locList = {'1','2'};                % location: 1: left screen, 2: right screen

epochedEEG = cell(length(code.mode),length(contList),length(locList),numTrial);
epochedIndex = zeros(length(code.mode),length(contList),length(locList));
for event_it = 1:length(EEG.event)
    tokenNames = regexp(EEG.event(event_it).type,expression,'names');
    if ~isempty(tokenNames)
        stimIndex = find(strcmp(tokenNames.stim,stimList));
        contIndex = find(strcmp(tokenNames.cont,contList));
        locIndex = find(strcmp(tokenNames.loc,locList));
        epochedIndex(stimIndex,contIndex,locIndex) = epochedIndex(stimIndex,contIndex,locIndex) + 1;
        
        dataIndex = floor(EEG.event(event_it).latency) + floor(epochWindow(1)*EEG.srate) : ...
            floor(EEG.event(event_it).latency) + floor(epochWindow(2)*EEG.srate);
        epochedEEG{stimIndex,contIndex,locIndex,epochedIndex(stimIndex,contIndex,locIndex)} = EEG.data(:,dataIndex);
    end
end

%% ------------------------------------------------------------------------
%                       EEG Analysis - deconvolution
% ------------------------------------------------------------------------

% select stimulation pattern: 
%   1: fmc-VEP (independent) 2: fmc-VEP (inverted) 3: fmc-VEP (chrome) 4: mseq
STIM = 1;
USE_EDGE = 0;   % 0: use raw stim code seq, 1: use difference / edge of the code seq

% select contrast level: 
%   1: +-4 level.  2: +-8 level 3: +-64 level
CONT = 3;

% select EEG channel ID(s) for analysis (1~34)
CH = 12; %[12,18:20]; %[5,10:14,18:20]; % [5,9:15,17:21,24:26]; %[12,17:21,24:26]; 

% define parameters
PLOT_RESULT = 0;    % plot frequency response of data, stim seq, impulse response, deconvolve filter
REFRESH = 60;       % screen refresh rate (in Hz)
len_trial = length(epochedEEG{1,1,1,1}(1,:));   % trial length (in samples)
windowing = hann(len_trial)';   % window applied for computing fft. Alternative option: ones(1,len_trial)

% leave-1-trial-out cross validation
acc = zeros(2,2);
rho = zeros(2,2);

for test_id = 1:numTrial  % define test trial ID 
    
    % store results for cross validation
    Gjw = cell(1,2);
    test_data = cell(1,2);
    S_true = cell(1,2);
    
    for LOC = 1:2   % for each location
        
        % extract all trials for the specified stim. pattern, contrast, location  
        data = [];
        for trial_id = 1:numTrial
            % average across all channels
            % [TODO] implement channel-wise deconvolution filter
            data = [data; mean(epochedEEG{STIM,CONT,LOC,trial_id}(CH,:),1)];
        end
        
        % extract training data (all but the test trial)
        X_train = data(setdiff(1:10,test_id),:);
        X_train = repmat(windowing,size(X_train,1),1) .* X_train;
        X_train_jw = fft(X_train,[],2);
               
        % extract raw code sequence
        if STIM == 1 || STIM == 3
            stim_seq = code.code_fmc{LOC};
        elseif STIM == 2
            stim_seq = code.code_fmc{LOC};
            if LOC == 2, stim_seq = ~stim_seq;  end
        elseif STIM == 4
            stim_seq = code.code_mseq{LOC};
        end
        
        % prepare stimulation code sequence
        S = zeros(1,len_trial);     % zero-pad at the end
        
        if ~USE_EDGE    % option 1: raw code sequence
            for it = 1:length(stim_seq)-1
                ind = floor( (it-1) * EEG.srate / REFRESH) + 1;
                S(ind) = 2*stim_seq(it)-1; % change to {-1, 1} to make it zero-mean
            end
        else    % option 2: difference / edge of code sequence
            edge_seq = diff(stim_seq);
            ind = floor((0:(length(edge_seq)-1)) * EEG.srate / REFRESH) + 1;
            S(ind) = edge_seq;
        end
        
        % compute fft of stimulation code sequence
        S = windowing .* S;
        S_jw = fft(S);         
 
        % compute impulse response: H = sum(XS*) / sum(SS* + sigma^2)
        regularized = 10;
        H_est_jw = sum(X_train_jw.*conj(S_jw)) ./ (size(X_train_jw,1)*(S_jw.*conj(S_jw)+regularized^2));
        H_est = ifft(H_est_jw);
        
        % compute deconvolution filter: G = 1 / H * ( |H|^2 / (|H|^2 + sigma^2) )
        regularized_deconv = 10;
        G_est_jw = 1./H_est_jw .* abs(H_est_jw).^2 ./ (abs(H_est_jw).^2 + regularized_deconv^2 );
        G_est = ifft(G_est_jw);
        
        % extract test data
        X_test = data(test_id,:);
        X_test = windowing .* X_test;
        X_test_jw = fft(X_test,[],2);

        % deconvolve signals
        S_est_jw = G_est_jw .* X_test_jw;
        S_est = ifft(S_est_jw);
                
        % plot result
        if PLOT_RESULT
            % average training trials - for plotting only
            X_avg = mean(X_train,1);
            X_avg = windowing .* X_avg;
            X_avg_jw = fft(X_avg);
            
            t = epochWindow(1):1/EEG.srate:epochWindow(2)+1/EEG.srate;
            f = 0:EEG.srate/len_trial:(EEG.srate-EEG.srate/len_trial);
            figure,
            subplot(2,2,1); plot(t,S); ylabel('s(t)'); xlabel('Time (sec)'); title('Stimulus sequence'); xlim([epochWindow(1) epochWindow(2)]);
            subplot(2,2,2); plot(f,abs(S_jw)); xlabel('Frequency (Hz)'); ylabel('|S(jw)|'); xlim([0 60]);
            subplot(2,2,3); plot(t,X_avg); ylabel('x(t)'); xlabel('Time (sec)'); title('EEG Data'); xlim([epochWindow(1) epochWindow(2)]);
            subplot(2,2,4); plot(f,abs(X_avg_jw)); xlabel('Frequency (Hz)'); ylabel('|X(jw)|'); xlim([0 60]);
            figure,
            subplot(2,2,1); plot(t,H_est); xlabel('Time (sec)'); ylabel('h(t)'); xlim([epochWindow(1) epochWindow(2)]);
            subplot(2,2,2); plot(f,abs(H_est_jw)); xlabel('Frequency (Hz)'); ylabel('|H(jw)|'); title('Estimated Impulse Resp.'); xlim([0 60]);
            subplot(2,2,3); plot(t,G_est); xlabel('Time (sec)'); ylabel('g(t)'); xlim([epochWindow(1) epochWindow(2)]);
            subplot(2,2,4); plot(f,abs(G_est_jw)); xlabel('Frequency (Hz)'); ylabel('|G(jw)|'); title('Estimated Deconv. Filter'); xlim([0 60]);
            figure,
            subplot(2,2,1); plot(t,S_est,'b');xlabel('Time (sec)'); ylabel('s(t)'); xlim([epochWindow(1) epochWindow(2)]);
            subplot(2,2,2); plot(f,abs(S_est_jw),'b'); xlabel('Frequency (Hz)'); ylabel('|S(jw)|'); title('Est. Stimulus Seq.'); xlim([0 60]);
            subplot(2,2,3); plot(t,S,'r');xlabel('Time (sec)'); ylabel('s(t)'); xlim([epochWindow(1) epochWindow(2)]);
            subplot(2,2,4); plot(f,abs(S_jw),'r'); xlabel('Frequency (Hz)'); ylabel('|S(jw)|'); title('True Stimulus Seq.'); xlim([0 60]);
        end
        
        Gjw{LOC} = G_est_jw;    % frequency response of deconvolution filter
        test_data{LOC} = X_test;% test trial
        S_true{LOC} = S;        % ground-true code sequence
       
    end
        
    % apply convolution filter to test data to obtain estimated code seq
    X_test_jw_1 = fft(test_data{1});
    X_test_jw_2 = fft(test_data{2});
    S_est_11 = ifft(Gjw{1} .* X_test_jw_1);
    S_est_12 = ifft(Gjw{1} .* X_test_jw_2);
    S_est_21 = ifft(Gjw{2} .* X_test_jw_1);
    S_est_22 = ifft(Gjw{2} .* X_test_jw_2);
    
    % compute cross-correlation
    rho11 = corrcoef(S_est_11,S_true{1});
    rho21 = corrcoef(S_est_21,S_true{1});
    rho12 = corrcoef(S_est_12,S_true{2});
    rho22 = corrcoef(S_est_22,S_true{2});
    rho(1,1) = rho(1,1) + rho11(1,2);
    rho(1,2) = rho(1,2) + rho12(1,2);
    rho(2,1) = rho(2,1) + rho21(1,2);
    rho(2,2) = rho(2,2) + rho22(1,2);
    
    % classify based on correlation
    if rho11(1,2) >= rho21(1,2)
        acc(1,1) = acc(1,1) + 1;
    else
        acc(1,2) = acc(1,2) + 1;
    end
    
    if rho22(1,2) >= rho12(1,2)
        acc(2,2) = acc(2,2) + 1;
    else
        acc(2,1) = acc(2,1) + 1;
    end
    
end
disp(acc)

%% ------------------------------------------------------------------------ 
%               TRCA-based VEP detection algorithm
% ------------------------------------------------------------------------- 

fprintf('Results of the ensemble TRCA-based method.\n');

% select stimulation pattern: 
%   1: fmc-VEP (independent) 2: fmc-VEP (inverted) 3: fmc-VEP (chrome) 4: mseq
STIM = 1;

% define parameters
numTrial = 10;      
is_ensemble = 1;    % use ensemble classifier (one classifier for each code seq)
channel_list = {[12,19], [12,18:20], [5,10:14,18:20], [5,9:15,17:21,24:26], [1:15,17:21,24:26], [1:32]};
data_len = 1.5;
data_offset = 0;

% define conditions for testing
cond_list = [1,2,3];        % contrast levels
var_list = channel_list;    % different channel combinations

numCH = zeros(1,length(channel_list));
mean_acc = zeros(length(cond_list),length(var_list));
mu_ci = zeros(length(cond_list),length(var_list),2);
for cond_i = 1:length(cond_list)    % constrast levels
    
    CONT = cond_list(cond_i);
    
    for var_i = 1:length(var_list)  % channels combinations
        
        CH = channel_list{var_i};
        numCH(var_i) = length(CH);
        
        data_range = data_offset + (1:1:floor(data_len*EEG.srate));
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
            %         fprintf('Trial %d: Accuracy = %2.2f%%\n', cv_i, accs(cv_i));
            
        end % loocv_i
        
        % Summarize results
        alpha_ci = 0.05;    % alpha for computing confidence interval
        ci = 1-alpha_ci;    % confidence interval
        [mu, ~, muci, ~] = normfit(accs, alpha_ci);
        fprintf('Mean accuracy = %2.2f %% (%2d%% CI: %2.2f - %2.2f %%)\n',...
            mu, ci, muci(1), muci(2));
        
        mean_acc(cond_i,var_i) = mu;
        mu_ci(cond_i,var_i,:) = muci;
    end
end

% plot classification results vs. number of channels for each contrast
figure, plot(numCH, mean_acc,'linewidth',2);
xlabel('Number of Channels'); ylabel('Cross validation accuracy (''%)');
set(gca,'fontsize',12)
ylim([0 100]); 
legend('Cont 4','Cont 8','Cont 64');


%% TRCA - # training trials
STIM = 4;
CONT = 1;
numTrial = 10;     
numCV = 20;         % cross validation trials
is_ensemble = 1;    % use ensemble classifier (one classifier for each code seq)

CH = [1:32]; %[5,9:15,17:21,24:26]; %[5,10:14,18:20]; %;
data_len = 3;
data_offset = 0;

nTrainTrial = 1:numTrial-1;

mean_acc = zeros(1,length(nTrainTrial));
mu_ci = zeros(2,length(nTrainTrial));

for numTrain = 1:length(nTrainTrial)
    
    data_range = data_offset + (1:1:floor(data_len*EEG.srate));
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
    accs = zeros(1,numCV);
    for cv_i = 1:numCV     % 10 times cross validation (selection with replacement)
        
        testID = randperm(10,10-numTrain);
        
        % Training stage
        traindata = trial_eeg;
        traindata(:, :, :, testID) = [];
        model = train_trca(traindata);
        
        % Test stage
        for test_i = 1:length(testID)
            testdata = squeeze(trial_eeg(:, :, :, testID(test_i)));
            estimated = test_trca(testdata, model, is_ensemble);
            
            % Evaluation
            is_correct = (estimated==labels);
            accs(cv_i) = accs(cv_i) + sum(is_correct);
            %         fprintf('Trial %d: Accuracy = %2.2f%%\n', cv_i, accs(cv_i));
        end
        accs(cv_i) = (accs(cv_i) / length(testID) / 2) * 100;
        
    end
    
    % Summarize results
    avg_acc = mean(accs);
    alpha_ci = 0.05;    % alpha for computing confidence interval
    ci = 1-alpha_ci;    % confidence interval
    [mu, ~, muci, ~] = normfit(accs, alpha_ci);
    fprintf('Mean accuracy = %2.2f %% (%2d%% CI: %2.2f - %2.2f %%)\n',...
        mu, ci, muci(1), muci(2));
    
    mean_acc(numTrain) = mu;
    mu_ci(:,numTrain) = muci;
    
end

% plot classification results vs. number of channels for each contrast
figure, plot(nTrainTrial, mean_acc,'linewidth',2);
xlabel('Number of Training Trials'); ylabel('Cross validation accuracy (''%)');
set(gca,'fontsize',12)
ylim([0 100]); 


%% Plot the task-related component from TRCA
STIM = 1;
CONT = 3;
CH = [5,10:14,18:20];
numTrial = 10;
is_ensemble = 1;
data_len = 3.2;
data_offset = 0; % linspace(0,floor((epochWindow(2)-epochWindow(1)-data_len)*EEG.srate), 25);

data_range = data_offset + (1:1:floor(data_len*EEG.srate));
NSAMP = length(data_range);

trial_eeg = zeros(2,length(CH),NSAMP,numTrial);
for loc_i = 1:2
    for tr_i = 1:10
        trial_eeg(loc_i,:,:,tr_i) = epochedEEG{STIM,CONT,loc_i,tr_i}(CH,data_range);
    end
end

% Leave-one-trial-out cross validation classification accuracy
labels = [1, 2];

% Training stage
traindata = trial_eeg;
model = train_trca(traindata);

source = cell(1,2);
for it = 1:model.num_targs
    source{it} = model.W(it,:) * squeeze(model.trains(it,:,:));
end

timeIndex = (1:length(source{1}))/EEG.srate;
figure,
subplot(3,1,1:2); hold on, 
plot(timeIndex, source{1}, 'b'); plot(timeIndex, source{2},'r'); legend('left','right')
subplot(3,1,3); plot(timeIndex, source{2}-source{1},'k')
xlabel('Time (sec)'); ylabel('Amplitude'); title('Source activity difference (R-L)')


%% ------------------------------------------------------------------------
%                 Ridge Regression for Decoding cVEP
% -------------------------------------------------------------------------
STIM = 1; 
CONT = 3;
CH = [12,19]; % [5,10:14,18:20]; % [12,18:20]; % [5,9:15,17:21,24:26]; %[12,17:21,24:26]; 
REFRESH = 60;

% define time window of EEG response to each code or edge
segWindow = round([-150, 50] * EEG.srate / 1000);    % samples

% extract timing of code sequence or that of rising / falling edges
on_fmc = cell(1,2);
off_fmc = cell(1,2);
rise_fmc = cell(1,2);
fall_fmc = cell(1,2);
flat_fmc = cell(1,2);
for loc_it = 1:2
    on_fmc{loc_it}  = find(code.code_fmc{loc_it} == 1);
    off_fmc{loc_it} = find(code.code_fmc{loc_it} == 0);
    rise_fmc{loc_it} = find((code.code_fmc{loc_it} - [0 code.code_fmc{loc_it}(1:end-1)]) == 1);
    fall_fmc{loc_it} = find((code.code_fmc{loc_it} - [0 code.code_fmc{loc_it}(1:end-1)]) == -1);
    flat_fmc{loc_it} = find((code.code_fmc{loc_it} - [0 code.code_fmc{loc_it}(1:end-1)]) == 0);
end

% extract EEG epochs time-locked to corresponding code sequence
X_on = [];
Y_on = [];
X_off = [];
Y_off = [];
for trial_id = 1:10
    for LOC = 1:2
        data = epochedEEG{STIM,CONT,LOC,trial_id}(CH,:);
        
        % extract EEG epochs time-locked to code 1
        for it = 1:length(on_fmc{loc_it})
            startIndex = round((on_fmc{loc_it}(it)-1) / REFRESH * EEG.srate) - segWindow(1);
            if startIndex+segWindow(2) <= length(data)
                temp_epoch = data(:, startIndex+1+segWindow(1) : startIndex+segWindow(2));
                X_on = [X_on; temp_epoch(:)' ];
                Y_on = [Y_on; 1];
            end
        end
        
        % extract EEG epochs time-locked to code 0
        for it = 1:length(off_fmc{loc_it})
            startIndex = round((off_fmc{loc_it}(it)-1) / REFRESH * EEG.srate) - segWindow(1);
            if startIndex+segWindow(2) <= length(data)
                temp_epoch = data(:, startIndex+1+segWindow(1) : startIndex+segWindow(2));
                X_off = [X_off; temp_epoch(:)' ];
                Y_off = [Y_off; -1];
            end
        end
        
    end
end

% separate training and testing data
CV_FOLD = 10;
alpha_list = 10.^(0:8);
[total_trial_on, feat_size] = size(X_on);
block_size = floor(total_trial / CV_FOLD);
acc_ls = zeros(1,CV_FOLD);
acc_ridge = zeros(length(alpha_list),CV_FOLD);
for cv_id = 1:CV_FOLD
    test_trial_id = (cv_id-1)*block_size+1 : cv_id*block_size;    
    train_trial_id = setdiff(1:total_trial, test_trial_id);
    X_test = [X_on(test_trial_id,:); X_off(test_trial_id,:)];
    Y_test = [Y_on(test_trial_id,:); Y_off(test_trial_id,:)];
    X_train = [X_on(train_trial_id,:); X_off(train_trial_id,:)];
    Y_train = [Y_on(train_trial_id,:); Y_off(train_trial_id,:)];
    
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


%% ------------------------------------------------------------------------
%                       Behavioral Results
% -------------------------------------------------------------------------

% load behavioral results
load('data\respMat_013119_pilot.mat')
code = respMat{1};
resp = respMat{2};

contrast_list = [4,8,64];

stat_percept = zeros(4,3,20);
stat_comfort = zeros(4,3,20);
index_mat = ones(4,3);
for it = 1:length(resp.stimuli)
   mode_index = 3 * resp.stimuli{it}(1) + resp.stimuli{it}(2) + 1;
   cont_index = find(resp.contrast(it)==contrast_list);
   
   stat_percept(mode_index,cont_index,index_mat(mode_index,cont_index)) = resp.rate_percept(it);
   stat_comfort(mode_index,cont_index,index_mat(mode_index,cont_index)) = resp.rate_comfort(it);
   index_mat(mode_index,cont_index) = index_mat(mode_index,cont_index) + 1;
end

mean_percept = mean(stat_percept,3);
std_percept = std(stat_percept,[],3);

mean_comfort = mean(stat_comfort,3);
std_comfort = std(stat_comfort,[],3);

figure
hold on
ha = bar(1:4,mean_percept);
% For each set of bars, find the centers of the bars, and write error bars
pause(0.1); %pause allows the figure to be created
for ib = 1:numel(ha)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = ha(ib).XData+ha(ib).XOffset;
    errorbar(xData,mean_percept(:,ib),std_percept(:,ib),'k.')
end
set(gca,'XTick',[1 2 3 4],'XTickLabel',{'Normal', 'Opposite', 'Chrome', 'Mseq'},'FontSize',12);
xlabel('Stimuli Types'); ylabel('Perceived Intensity Rating');

figure
hold on
hb = bar(1:4,mean_comfort);
% For each set of bars, find the centers of the bars, and write error bars
pause(0.1); %pause allows the figure to be created
for ib = 1:numel(hb)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hb(ib).XData+hb(ib).XOffset;
    errorbar(xData,mean_comfort(:,ib),std_comfort(:,ib),'k.')
end
ylim([2 6]);
set(gca,'XTick',[1 2 3 4],'XTickLabel',{'Normal', 'Opposite', 'Chrome', 'Mseq'},'FontSize',12);
xlabel('Stimuli Types'); ylabel('Comfortability Rating');

%% ------------------------------------------------------------------------
%       Extract template of impulse response to rising / falling edges
% -------------------------------------------------------------------------

STIM = 1; 
CONT = 3;
CH = [5,10:14,18:20]; % [12,18:20]; % [5,9:15,17:21,24:26]; %[12,17:21,24:26]; 
REFRESH = 60;

% extract and visualize the impulse response to rising edges
% lenSegment = 0.5;   % sec
segWindow = round([-200, 400] * EEG.srate / 1000);    % samples
data_seg_rise = []; % zeros(numSegment, round(lenSegment*EEG.srate));
data_seg_fall = []; % zeros(numSegment, round(lenSegment*EEG.srate));
data_seg_flat = [];

% extract code sequence and rising / falling edges
rise_fmc = cell(1,2);
fall_fmc = cell(1,2);
flat_fmc = cell(1,2);
for loc_it = 1:2
    rise_fmc{loc_it} = find((code.code_fmc{loc_it} - [0 code.code_fmc{loc_it}(1:end-1)]) == 1);
    fall_fmc{loc_it} = find((code.code_fmc{loc_it} - [0 code.code_fmc{loc_it}(1:end-1)]) == -1);
    flat_fmc{loc_it} = find((code.code_fmc{loc_it} - [0 code.code_fmc{loc_it}(1:end-1)]) == 0);
end


for trial_id = 1:10
    for LOC = 1:2
        data = mean(epochedEEG{STIM,CONT,LOC,trial_id}(CH,:),1);
        
        % extract the impulse response to rising edges
        for it = 1:length(rise_fmc{LOC})
            startIndex = round((rise_fmc{LOC}(it)-1) / REFRESH * EEG.srate) - segWindow(1);
            if startIndex+segWindow(2) <= length(data)
                data_seg_rise = [data_seg_rise; ...
                    data( startIndex+1+segWindow(1) : startIndex+segWindow(2)) ];
            end
        end
    
        % extract the impulse response to falling edges
        for it = 1:length(fall_fmc{LOC})
            startIndex = round((fall_fmc{LOC}(it)-1) / REFRESH * EEG.srate) - segWindow(1);
            if startIndex+segWindow(2) <= length(data)
                data_seg_fall = [data_seg_fall; ...
                    data( startIndex+1+segWindow(1) : startIndex+segWindow(2)) ];
            end
        end
        
        % extract the impulse response to no edges
        for it = 1:length(flat_fmc{LOC})
            startIndex = round((flat_fmc{LOC}(it)-1) / REFRESH * EEG.srate) - segWindow(1);
            if startIndex+segWindow(2) <= length(data)
                data_seg_flat = [data_seg_flat; ...
                    data( startIndex+1+segWindow(1) : startIndex+segWindow(2)) ];
            end
        end
    end
end
avg_data_seg_rise = mean(data_seg_rise,1);
avg_data_seg_fall = mean(data_seg_fall,1);
avg_data_seg_flat = mean(data_seg_flat,1);
avg_data_seg = mean([data_seg_rise; -data_seg_fall],1);

% plot the average impulse response
figure, subplot(2,1,1); hold on, 
plot((segWindow(1)+1:segWindow(2))./EEG.srate*1000,avg_data_seg_rise,'b')
plot((segWindow(1)+1:segWindow(2))./EEG.srate*1000,avg_data_seg_fall,'r')
xlim([-200 400]); legend('rising edge','falling edge'); title('Avg Impulse Response (IR)')
subplot(2,1,2);
plot((segWindow(1)+1:segWindow(2))./EEG.srate*1000, avg_data_seg, 'k')
xlabel('Time (msec)'); ylabel('Amplitude (\muV)'); xlim([-200 400]); title('IR(rising)-IR(falling)')
% plot((segWindow(1)+1:segWindow(2))./EEG.srate*1000,avg_data_seg_flat,'k')

% plot the power spectra 
[prise,frise] = periodogram(avg_data_seg_rise,[],[],EEG.srate);
[pfall,ffall] = periodogram(avg_data_seg_fall,[],[],EEG.srate);
[pdiff,fdiff] = periodogram(avg_data_seg,[],[],EEG.srate);

figure, subplot(2,1,1); hold on, 
plot(frise,20*log10(abs(prise)),'b')
plot(ffall,20*log10(abs(pfall)),'r')
legend('rising edge','falling edge'); title('Avg Impulse Response (IR)')
subplot(2,1,2);
plot(fdiff,20*log10(abs(pdiff)),'k')
xlabel('Frequency (Hz)'); ylabel('Spectral Power (dB)'); title('IR(rising)-IR(falling)')


%% Template matching correlation over time
testCH = [5,10:14,18:20]; % [12,18:20]; 
LOC = 1;

% extract and visualize the impulse response to rising edges
template_rising = avg_data_seg_rise;
template_falling = avg_data_seg_fall;
template_flat = avg_data_seg_flat;
template_diff = avg_data_seg;

len_template = length(avg_data_seg_rise);
len_data = size(epochedEEG{1,1,1,1},2);
corr_rising = zeros(1,len_data-len_template+1);
corr_falling= zeros(1,len_data-len_template+1);
corr_flat = zeros(1,len_data-len_template+1);
corr_diff = zeros(1,len_data-len_template+1);

data = zeros(10,len_data);
for trial_id = 1:10
    data(trial_id,:,:) = mean(epochedEEG{STIM,CONT,LOC,trial_id}(testCH,:),1);
end
data = mean(data,1);
  
for it = 1:len_data-len_template+1
    corr_rising(it) = corr(data(it:it+len_template-1)',template_rising');
    corr_falling(it) = corr(data(it:it+len_template-1)',template_falling');
    corr_flat(it) = corr(data(it:it+len_template-1)',template_flat');
end

% plot correlation of template matching over time
figure, hold on,
plot(corr_rising,'b');
plot(corr_falling,'r');
plot(corr_flat,'k');
xlabel('Samples'); ylabel('Correlation'); legend('Rising','Falling','Flat'); set(gca,'fontsize',12);

numCode = floor((len_data-len_template+1)/EEG.srate*REFRESH);
maxCorr = zeros(floor(EEG.srate/REFRESH),numCode-1);
decodedSeq = zeros(floor(EEG.srate/REFRESH),numCode-1);
for shift = 1:floor(EEG.srate/REFRESH)
    for code_id = 1:numCode-1
        codeIndex = round((code_id-1)*EEG.srate/REFRESH) + shift + (0:floor(EEG.srate/REFRESH-1));
        [maxCorr(shift,code_id),decodedSeq(shift,code_id)] = max( ...
            mean([corr_rising(codeIndex);corr_flat(codeIndex);corr_falling(codeIndex)],2));
    end
end
[~,opt_shift] = max(median(maxCorr,2));
codeMatching = decodedSeq(opt_shift,:);
codeConfidence = maxCorr(opt_shift,:);
% figure, plot(codeMatching);

codeSeq = zeros(1,length(codeMatching));
for it = 1:length(codeMatching)
    if codeMatching(it) == 1  % rising edge
        codeSeq(it) = 1;
    elseif codeMatching(it) == 3  % falling edge
        codeSeq(it) = 0;
    elseif codeMatching(it) == 2 && it > 1 % code: 2 no change
        codeSeq(it) = codeSeq(it-1);
    else
        codeSeq(it) = 0;
    end
end

figure, 
subplot(2,1,1); imagesc(codeSeq);
subplot(2,1,2); imagesc(codeConfidence);

%% Misc.
LOC = 1;
errorRate = zeros(1,length(code.code_fmc{LOC})-length(decodedSeq));
for it = 1:length(code.code_fmc{LOC})-length(decodedSeq)
    errorRate(it) = mean(code.code_fmc{LOC}(it:it+length(decodedSeq)-1) ~= codeSeq);
end

[minErrorRate,opt_align] = min(errorRate);

figure,
subplot(3,1,1); imagesc(code.code_fmc{LOC}(opt_align:opt_align+length(decodedSeq)-1)); title('Actual Code');
subplot(3,1,2); imagesc(codeSeq); title('Estimated Code')
subplot(3,1,3); imagesc(code.code_fmc{LOC}(opt_align:opt_align+length(decodedSeq)-1)~=codeSeq); title('Error');
xlabel('Time (samples)'); set(gca,'fontsize',12)

%{
[peaks_rising,locs_rising] = findpeaks(corr_rising,'MinPeakHeight',0,'MinPeakDistance',8);
[peaks_falling,locs_falling] = findpeaks(corr_falling,'MinPeakHeight',0,'MinPeakDistance',8);
[peaks_flat,locs_flat] = findpeaks(corr_flat,'MinPeakHeight',0,'MinPeakDistance',8);

figure, hold on,
stem(locs_rising,peaks_rising,'b');
stem(locs_falling,peaks_falling,'r');
stem(locs_flat,peaks_flat,'k');
legend('Rising','Falling','Flat'); xlabel('Samples'); ylabel('Correlation'); set(gca,'fontsize',12);

[sortedLoc, sortedIndex] = sort([locs_rising, locs_falling]);
decoder = [ones(1,length(locs_rising)),zeros(1,length(locs_falling))];
decoder = decoder(sortedIndex);
figure, stem(sortedLoc,decoder);

tmp = mod(sortedLoc, EEG.srate/REFRESH);
figure, hist(tmp,20);
%}


%% Generate data matrix
len_trial = length(epochedEEG{1,1,1,1}(1,:));   % trial length (in samples)
SRATE = 512;
REFRESH = 60;
OFFSET = 50; % msec

stim_data = cell(size(epochedEEG));
for STIM = 1:4
    for LOC = 1:2
        % extract raw code sequence
        if STIM == 1 || STIM == 3
            code = stim_seq.fmc{LOC};
        elseif STIM == 2
            code = stim_seq.fmc_inv{LOC};
        elseif STIM == 4
            code = stim_seq.mseq{LOC};
        end
        
        % prepare stimulation code sequence
        S = zeros(1,len_trial);     % zero-pad at the end
        S_edge = zeros(1,len_trial);     % zero-pad at the end
        ini_index = floor(OFFSET / 1000 * SRATE);
        for it = 1:length(code)-1
            ind = floor( (it-1) * SRATE / REFRESH) + 1 + ini_index;
            S(ind) = 2*code(it)-1; % change to {-1, 1} to make it zero-mean
        end
        edge_seq = diff(code);
        ind = floor((0:(length(edge_seq)-1)) * SRATE / REFRESH) + 1 + ini_index;
        S_edge(ind) = edge_seq;
        
        for CONT = 1:3
            for TRIAL = 1:10
                stim_data{STIM,CONT,LOC,TRIAL} = [S; S_edge];
            end
        end
    end
end

save('stim_data.mat','stim_data')

