
REFRESH = 60;   % Hz
SRATE = 512;    % Hz
T = 3;  % data length (sec)
N = T*REFRESH;   % stimulus length

% generate stimulus sequence
rand_seq = rand(1,N) > 0.5;
stim_seq = diff(rand_seq);
S = zeros(1,T*SRATE);
ind = floor((0:1:(N-2)) * SRATE / REFRESH) + 1;
S(ind) = stim_seq;
S_jw = fft(S);

% generate impulse response
tao = 10;
freq = 10;
t = 0:1/SRATE:(T-1/SRATE); 
H = exp(-tao*t) .* sin(2*pi*freq*t);
H = H + awgn(H,20);
H_jw = fft(H);

% simulated signal
X = conv(H, S, 'same');
% X = X(1:length(t));
X_jw = fft(X);

% plot simulation results
f = 0:1/T:(SRATE-1/T);
figure, 
subplot(3,2,1); plot(t,S); ylabel('s(t)'); xlabel('Time (sec)'); title('stimulus sequence');
subplot(3,2,2); plot(f,abs(S_jw)); xlabel('Frequency (Hz)'); ylabel('|S(jw)|');
subplot(3,2,3); plot(t,H); ylabel('h(t)'); xlabel('Time (sec)'); title('impulse response');
subplot(3,2,4); plot(f,abs(H_jw)); xlabel('Frequency (Hz)'); ylabel('|H(jw)|');
subplot(3,2,5); plot(t,X); ylabel('x(t)'); xlabel('Time (sec)'); title('simulated data');
subplot(3,2,6); plot(f,abs(X_jw)); xlabel('Frequency (Hz)'); ylabel('|X(jw)|');

%% 
% option 1: H = X / S (same as option 2)
H_jw_est_1 = X_jw ./ S_jw;                % zero at 0 Hz
H_jw_est_1(1) = 0;                        % set zero frequency to 0
H_est_1 = ifft(H_jw_est_1);

% option 2: H = XX* / SX*
H_jw_est = (X_jw .* conj(X_jw)) ./ (S_jw .* conj(X_jw));
H_jw_est(1) = 0;
H_est = ifft(H_jw_est);

% plot results
figure,
subplot(2,2,1); plot(f,abs(H_jw_est),'b'); xlabel('Frequency (Hz)'); ylabel('|H(jw)|'); title('Estimated Impulse Resp.')
subplot(2,2,2); plot(f,abs(H_jw),'r'); xlabel('Frequency (Hz)'); ylabel('|H(jw)|'); title('True Impulse Resp.')
subplot(2,2,3); plot(t,H_est,'b'); xlabel('Time (sec)'); ylabel('Est.');
subplot(2,2,4); plot(t,H,'r');xlabel('Time (sec)'); ylabel('True');


%%
% compute deconvolve filter
G_jw_est = 1./H_jw_est;
G_jw_est(1) = 0;
G_est = ifft(G_jw_est);

% deconvlve signals
S_jw_est = G_jw_est .* X_jw;
S_est = ifft(S_jw_est);

% plot results
figure,
subplot(3,2,1); plot(f,abs(G_jw_est),'k'); xlabel('Frequency (Hz)'); ylabel('|G(jw)|'); title('Deconv. Filter')
subplot(3,2,2); plot(t,G_est,'k'); xlabel('Time (sec)'); ylabel('g(t)');
subplot(3,2,3); plot(f,abs(S_jw_est),'b'); xlabel('Frequency (Hz)'); ylabel('|S(jw)|'); title('Est. Stimulus Seq.')
subplot(3,2,4); plot(t,S_est,'b');xlabel('Time (sec)'); ylabel('s(t)'); ylim([-1 1]);
subplot(3,2,5); plot(f,abs(S_jw),'r'); xlabel('Frequency (Hz)'); ylabel('|S(jw)|'); title('True Stimulus Seq.')
subplot(3,2,6); plot(t,S,'r');xlabel('Time (sec)'); ylabel('s(t)');


%% 
% apply hanning window
X_win = hann(length(X))' .* X;
X_win_jw = fft(X_win);

% deconvolution
[Q,R] = deconv(X,S);
