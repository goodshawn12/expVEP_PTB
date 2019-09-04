
STIM = 3;
CONT = 3;
CH = 12; 

figure, hold on,
for LOC = 1:2
for trial_id = 1:10

% % extract all trials for the specified stim. pattern, contrast, location
% data = [];
% num_trial= 10;
% for trial_id = 1:num_trial
%     % average across all channels and stack trials as rows in data
%     % [TODO] implement channel-wise deconvolution filter
%     data = [data; mean(epochedEEG{STIM,CONT,LOC,trial_id}(:,:),1)];
% end
% sample_data = mean(data);

sample_data = mean(epochedEEG{STIM,CONT,LOC,trial_id}(CH,:),1);

% figure, plot(sample_data)
% psd_data = fftshift(fft(sample_data));
% f = linspace(-256, 256, length(sample_data));
% figure, plot(f, abs(psd_data))

srate = 512;
window_len = 1;
window = hann(window_len*srate);
noverlap = floor(length(window)/2);
nfft = 2.^ceil(log2(length(window)));
[Pxx,F] = pwelch(sample_data,window,noverlap,nfft,srate);
subplot(10,2,(trial_id-1)*2+LOC), plot(F,Pxx,'linewidth',2); xlim([5 60]);

% srate = 512;
% window_len = 0.5;
% window = hann(window_len*srate);
% noverlap = floor(length(window)*0.9);
% nfft = 2.^ceil(log2(length(window)));
% [s,f,t,p] = spectrogram(sample_data,window,noverlap,nfft,srate);
% figure, surf(t,f,db(p))
% ylim([0 60])
end
end
