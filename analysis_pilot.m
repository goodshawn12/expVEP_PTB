
% addpath EEGLAB
cd('C:\Users\shawn\Desktop\NSF_BIC\expVEP_PTB')


%% ------------------------------------------------------------------------
%     EEG Analysis - Average Impulse Response to Rising / Falling Edges
% -------------------------------------------------------------------------
filename = 'expVEP_013119_EEG.set';
filepath = 'C:\Users\shawn\Desktop\NSF_BIC\expVEP_PTB\data\';
EEG = pop_loadset([filepath filename]);

load('data\respMat_013119_pilot.mat')
code = respMat{1};

% extract epoch
expression = 'Stim(?<stim>\d+)_CT(?<cont>\d+)_Loc(?<loc>\d+)';
stimList = {'00','01','02','10'};
contList = {'4','8','64'};
locList = {'1','2'};
epochedEEG = cell(length(code.mode),3,2,10);
epochedIndex = zeros(length(code.mode),3,2);
epochWindow = [-0.2, 3];   % sec
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

%% Deconvolution Analysis
STIM = 1; 
CONT = 3;
CH = 12; %[12,18:20]; %[5,10:14,18:20]; % [5,9:15,17:21,24:26]; %[12,17:21,24:26]; 
REFRESH = 60;

acc = zeros(2,2);
rho = zeros(2,2);
for test_id = 1:10
    gt = cell(1,2);
    test_data = cell(1,2);
    S_true = cell(1,2);
    for LOC = 1:2
        data = [];
        for trial_id = 1:10
            data = [data; mean(epochedEEG{STIM,CONT,LOC,trial_id}(CH,:),1)];
        end
        
        X_train = data(setdiff(1:10,test_id),:);
        train_jw = fft(X_train,[],2);
        
        X_test = data(test_id,:);
        X = mean(X_train,1);
        window = hann(length(X))'; % ones(1,length(X)); % 
        X = window .* X;
        X_jw = fft(X);
        
        % extract stored code
        if STIM == 1 || STIM == 3
            stim_seq = code.code_fmc{LOC};
        elseif STIM == 2
            stim_seq = code.code_fmc{LOC};
            if LOC == 2, stim_seq = ~stim_seq;  end
        elseif STIM == 4
            stim_seq = code.code_mseq{LOC};
        end
        S = zeros(1,length(X));     % zero-pad at the end
        
        % option 1: raw code
%         for it = 1:length(stim_seq)-1
%             ind = floor( (it-1) * EEG.srate / REFRESH) + 1;
%             S(ind) = 2*stim_seq(it)-1;        % change to {-1, 1} to make it zero-mean
%         end
        
        % option 2: difference code
        edge_seq = diff(stim_seq);
        ind = floor((0:(length(edge_seq)-1)) * EEG.srate / REFRESH) + 1;
        S(ind) = edge_seq;
        
        %
        S = window .* S;
        S_jw = fft(S); 
        
        % compute impulse response
        % option A: H = X / S
%         H_jw_est = X_jw ./ S_jw;                % zero at 0 Hz
%         H_jw_est(1) = 0;
%         H_est = ifft(H_jw_est);
 
        % option B: H = sum(X X*) / sum(S X*)
        noise = 50;
        H_jw_est = sum(train_jw.*conj(S_jw)) ./ (size(train_jw,1)*(S_jw.*conj(S_jw)+noise^2));
        H_est = ifft(H_jw_est);
       
        % compute deconvolve filter
        G_jw_est = 1./H_jw_est;
%         G_jw_est(1) = 0;
        G_est = ifft(G_jw_est);
        
        % deconvlve signals
        X_jw_test = fft(X_test);
        S_jw_est = G_jw_est .* X_jw_test;
        S_est = ifft(S_jw_est);
                
        % plot result
        t = epochWindow(1):1/EEG.srate:epochWindow(2)+1/EEG.srate;
        f = 0:EEG.srate/length(X):(EEG.srate-EEG.srate/length(X));
        figure,
        subplot(2,2,1); plot(t,S); ylabel('s(t)'); xlabel('Time (sec)'); title('Stimulus sequence'); xlim([epochWindow(1) epochWindow(2)]);
        subplot(2,2,2); plot(f,abs(S_jw)); xlabel('Frequency (Hz)'); ylabel('|S(jw)|'); xlim([0 60]);
        subplot(2,2,3); plot(t,X); ylabel('x(t)'); xlabel('Time (sec)'); title('EEG Data'); xlim([epochWindow(1) epochWindow(2)]);
        subplot(2,2,4); plot(f,abs(X_jw)); xlabel('Frequency (Hz)'); ylabel('|X(jw)|'); xlim([0 60]);
        figure,
        subplot(2,2,1); plot(t,H_est); xlabel('Time (sec)'); ylabel('h(t)'); xlim([epochWindow(1) epochWindow(2)]);
        subplot(2,2,2); plot(f,abs(H_jw_est)); xlabel('Frequency (Hz)'); ylabel('|H(jw)|'); title('Estimated Impulse Resp.'); xlim([0 60]);
        subplot(2,2,3); plot(t,G_est); xlabel('Time (sec)'); ylabel('g(t)'); xlim([epochWindow(1) epochWindow(2)]);
        subplot(2,2,4); plot(f,abs(G_jw_est)); xlabel('Frequency (Hz)'); ylabel('|G(jw)|'); title('Estimated Deconv. Filter'); xlim([0 60]);
        figure,
        subplot(2,2,1); plot(t,S_est,'b');xlabel('Time (sec)'); ylabel('s(t)'); xlim([epochWindow(1) epochWindow(2)]);
        subplot(2,2,2); plot(f,abs(S_jw_est),'b'); xlabel('Frequency (Hz)'); ylabel('|S(jw)|'); title('Est. Stimulus Seq.'); xlim([0 60]);
        subplot(2,2,3); plot(t,S,'r');xlabel('Time (sec)'); ylabel('s(t)'); xlim([epochWindow(1) epochWindow(2)]);
        subplot(2,2,4); plot(f,abs(S_jw),'r'); xlabel('Frequency (Hz)'); ylabel('|S(jw)|'); title('True Stimulus Seq.'); xlim([0 60]);
        
        gt{LOC} = G_est;
        test_data{LOC} = X_test;
        S_true{LOC} = S;
        
        % testing
        H_jw_trunc = H_jw_est;
        trunc_indx = find(f>60,1);
        H_jw_trunc(trunc_indx:end-trunc_indx) = 0;
        H_trunc_est = ifft(H_jw_trunc);
        
        figure, plot(f,abs(H_jw_est),'b',f,abs(H_jw_trunc),'r');
        figure, plot(t,H_est,'b',t,H_trunc_est,'r');

        regularize = 0.1;
        G_jw_reg = 1./ (H_jw_est + regularize^2);
        G_est_reg = ifft(G_jw_reg);
        
        figure, plot(f,abs(G_jw_est),'b',f,abs(G_jw_reg),'r');
        figure, plot(t,G_est,'b',t,G_est_reg,'r');
        
        S_jw_reg = G_jw_reg .* X_jw_test;
        S_est_reg = ifft(S_jw_reg);
        
        figure, plot(f,abs(S_jw_est),'b',f,abs(S_jw_reg),'r');
        figure, plot(t,S_est,'b',t,S_est_reg,'r');

    end
    
    % autocorrelation
    rho = corrcoef(S_est,S); disp(rho(1,2));
    rho = corrcoef(S_est_reg,S); disp(rho(1,2));
    

    
    % classification using deconvolution filter
    G_est_jw_1 = fft(gt{1});
    G_est_jw_2 = fft(gt{2});
    X_test_jw_1 = fft(test_data{1});
    X_test_jw_2 = fft(test_data{2});
    S_est_jw_11 = G_est_jw_1 .* X_test_jw_1;
    S_est_jw_12 = G_est_jw_1 .* X_test_jw_2;
    S_est_jw_21 = G_est_jw_2 .* X_test_jw_1;
    S_est_jw_22 = G_est_jw_2 .* X_test_jw_2;
    S_est_11 = ifft(S_est_jw_11);
    S_est_12 = ifft(S_est_jw_12);
    S_est_21 = ifft(S_est_jw_21);
    S_est_22 = ifft(S_est_jw_22);
    
    rho11 = corrcoef(S_est_11,S_true{1});
    rho21 = corrcoef(S_est_21,S_true{1});
    rho12 = corrcoef(S_est_12,S_true{2});
    rho22 = corrcoef(S_est_22,S_true{2});
    rho(1,1) = rho(1,1) + rho11(1,2);
    rho(1,2) = rho(1,2) + rho12(1,2);
    rho(2,1) = rho(2,1) + rho21(1,2);
    rho(2,2) = rho(2,2) + rho22(1,2);
    
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
    
    % S_est_conv_l1 = conv(gt{1},test_data{1});
    % S_est_conv_21 = conv(gt{2},test_data{1});
    % S_est_conv_l2 = conv(gt{1},test_data{2});
    % S_est_conv_22 = conv(gt{2},test_data{2});
end
disp(acc)

%% select channel and average over trials - template of impulse response
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

%%
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



%% ------------------------------------------------------------------------ 
%       Performing the TRCA-based VEP detection algorithm
% ------------------------------------------------------------------------- 

fprintf('Results of the ensemble TRCA-based method.\n');

% Preparing data
STIM = 1;
CONT = 3;
CH = [5,10:14,18:20];
NTRIAL = 10;
is_ensemble = 1;

channel_list = {[12,19], [12,18:20], [5,10:14,18:20], [5,9:15,17:21,24:26], [1:15,17:21,24:26], [1:32]};

data_len = 1.5;
data_offset = 0; % linspace(0,floor((epochWindow(2)-epochWindow(1)-data_len)*EEG.srate), 25);

cond_list = [1,2,3];
var_list = channel_list;
numCH = zeros(1,length(channel_list));

mean_acc = zeros(length(cond_list),length(var_list));
mu_ci = zeros(length(cond_list),length(var_list),2);
for cond_i = 1:length(cond_list)
    
    CONT = cond_list(cond_i);
    
    for var_i = 1:length(var_list)
        
        CH = channel_list{var_i};
        numCH(var_i) = length(CH);
        
        data_range = data_offset + (1:1:floor(data_len*EEG.srate));
        NSAMP = length(data_range);
        
        trial_eeg = zeros(2,length(CH),NSAMP,NTRIAL);
        for loc_i = 1:2
            for tr_i = 1:10
                trial_eeg(loc_i,:,:,tr_i) = epochedEEG{STIM,CONT,loc_i,tr_i}(CH,data_range);
            end
        end
        
        % Leave-one-trial-out cross validation classification accuracy
        labels = [1, 2];
        for cv_i = 1:1:NTRIAL
            
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
        
        % Summarize
        alpha_ci = 0.05;
        ci = 1-alpha_ci;
        [mu, ~, muci, ~] = normfit(accs, alpha_ci);
        fprintf('Mean accuracy = %2.2f %% (%2d%% CI: %2.2f - %2.2f %%)\n',...
            mu, ci, muci(1), muci(2));
        
        mean_acc(cond_i,var_i) = mu;
        mu_ci(cond_i,var_i,:) = muci;
    end
end
% figure, plot(data_offset_list/EEG.srate + epochWindow(1), mean_acc,'linewidth',2);
% xlabel('Window Delay (sec)'); ylabel('Cross validation accuracy (''%)');
% set(gca,'fontsize',12)

% figure, plot(data_len_list, mean_acc,'linewidth',2);
% xlabel('Trial length (sec)'); ylabel('Cross validation accuracy (''%)');
% set(gca,'fontsize',12)
% ylim([0 100]); 
% legend('FMC-idp','FMC-opp','MSEQ');
% legend('FMC-idp','FMC-opp','FMC-chm','MSEQ');
% legend('Cont 4','Cont 8','Cont 64');

figure, plot(numCH, mean_acc,'linewidth',2);
xlabel('Number of Channels'); ylabel('Cross validation accuracy (''%)');
set(gca,'fontsize',12)
ylim([0 100]); 
legend('Cont 4','Cont 8','Cont 64');

%% Plot the task-related component from TRCA
STIM = 1;
CONT = 3;
CH = [5,10:14,18:20];
NTRIAL = 10;
is_ensemble = 1;
data_len = 3.2;
data_offset = 0; % linspace(0,floor((epochWindow(2)-epochWindow(1)-data_len)*EEG.srate), 25);

data_range = data_offset + (1:1:floor(data_len*EEG.srate));
NSAMP = length(data_range);

trial_eeg = zeros(2,length(CH),NSAMP,NTRIAL);
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
%                   Behavioral Results
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