% ------------------------------------------------------------------------
%                       EEG Analysis - deconvolution
% ------------------------------------------------------------------------
function vep_deconv(epochedEEG,stim_data)

% select stimulation pattern:
%   1: fmc-VEP (independent) 2: fmc-VEP (inverted) 3: fmc-VEP (chrome) 4: mseq
STIM = 1;

% select contrast level:
%   1: +-4 level.  2: +-8 level 3: +-64 level
CONT = 3;

% select location:
%   1: left screen. 2: right screen
LOC = 1;

% select EEG channel ID(s) for analysis (1~34)
CH = 12; %[12,18:20]; %[5,10:14,18:20]; % [5,9:15,17:21,24:26]; %[12,17:21,24:26];

% define parameters
PLOT_RESULT = 1;    % plot frequency response of data, stim seq, impulse response, deconvolve filter
num_trial = 10;
len_trial = length(epochedEEG{1,1,1,1}(1,:));   % trial length (in samples)
SRATE = 512;
epochWindow = [-0.2 3]; % sec

% define window applied to data for computing fft.
windowing = hann(len_trial)';  % Alternative option: ones(1,len_trial)

% test trial index (1~10)
test_id = 10;


% extract all trials for the specified stim. pattern, contrast, location
data = [];
for trial_id = 1:num_trial
    % average across all channels and stack trials as rows in data
    % [TODO] implement channel-wise deconvolution filter
    data = [data; mean(epochedEEG{STIM,CONT,LOC,trial_id}(CH,:),1)];
end

% extract training data (all but the test trial)
X_train = data(setdiff(1:num_trial,test_id),:);
X_train = repmat(windowing,size(X_train,1),1) .* X_train;
X_train_jw = fft(X_train,[],2);

% extract stimulus sequence data
S_raw = stim_data{STIM,CONT,LOC,1}(1,:);
S_edge = stim_data{STIM,CONT,LOC,1}(2,:);
S = S_edge;

% compute fft of stimulation code sequence
S = windowing .* S;
S_jw = fft(S);

% compute impulse response: H = sum(XS*) / sum(SS* + sigma^2)
regularized = 10;   % avoid division of small value in high-freq domain
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
    
    t = epochWindow(1) : (1/SRATE) : (epochWindow(2)+1/SRATE);
    f = 0:SRATE/len_trial:(SRATE-SRATE/len_trial);
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
end
